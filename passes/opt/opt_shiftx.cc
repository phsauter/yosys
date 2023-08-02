/*
 *  yosys -- Yosys Open SYnthesis Suite
 *
 *  Copyright (C) 2012  Claire Xenia Wolf <claire@yosyshq.com>
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include "kernel/yosys.h"
#include "kernel/sigtools.h"
#include "kernel/ffinit.h"

USING_YOSYS_NAMESPACE
PRIVATE_NAMESPACE_BEGIN

struct OptShiftxWorker
{
	struct shiftxnode_t {
		Cell *shiftx;

		SigSpec select_idx;

		int block_size;
		int block_offset;
		bool shift_left;
	};

	RTLIL::Design *design;
	RTLIL::Module *module;
	SigMap sigmap;
	int min_shift;
	bool napot;
	bool verbose;

	dict<SigBit, Cell *> sigbit_to_driver; // TODO: I can't use SigChunk because it has no hash func, odd?
	std::vector<Cell *> shiftx_cells;
	std::vector<shiftxnode_t> candiadates;

	OptShiftxWorker(RTLIL::Design *design, RTLIL::Module *module, int min_shift, bool napot, bool verbose) :
			design(design), module(module), sigmap(module), min_shift(min_shift), napot(napot), verbose(verbose) { }


	void index()
	{
		// extract important info from module
		for (auto cell : module->cells())
		{
			// track drivers and connected wires of all cells
			for (auto &conn : cell->connections())
				if (cell->output(conn.first))
					for (auto bit : sigmap(conn.second))
						sigbit_to_driver[bit] = cell;


			// track selected shift/shiftx cells
			if(cell->type.in(ID($shift), ID($shiftx))) {
				if (design->selected(module, cell))
					shiftx_cells.push_back(cell);
			}
		}

		if(verbose) log("  %s: indexed, found %ld $shift/$shiftx cells\n", log_id(module), shiftx_cells.size());
	}

	void check_shiftx(Cell *shiftx)
	{
		shiftxnode_t shiftxnode;
		shiftxnode.shiftx = shiftx;

		SigSpec B = sigmap(shiftx->getPort(ID::B)); // shift-amount (selection)
		Cell *prev = nullptr;

		if(verbose) log("  %s: following cell: %s.\n", log_id(module), log_id(shiftx));

		// in principle at any stage the signal can be concatenated with any value, not just zero
		// so tracking the number of const bits is likely insufficent
		int lsb_const_bits = 0;
		for(auto bit : B.bits()) {
			if(bit.is_wire())
				prev = sigbit_to_driver.at(bit, nullptr);
			else
				lsb_const_bits++;

			if(prev != nullptr)
				break;
		}

		int block_size = 0;
		// used to keep track of the offset of the first block
		bool block_offset_neg = false;
		int block_offset = 0;

		// only allow constant block_offsets (+/-) before seing mul
		while(prev && prev->type.in(ID($add), ID($sub), ID($neg)))
		{
			if(verbose) log("  %s: in input cone %s\n", log_id(module), log_id(prev));

			// negating only has one input, try to follow it
			SigSpec prev_A = sigmap(prev->getPort(ID::A));
			if(prev->type == ID($neg)) {
				block_offset_neg = !block_offset_neg;
				prev = nullptr;
				for(auto bit : prev_A.bits()) {
					if(bit.is_wire())
						prev = sigbit_to_driver.at(bit, nullptr);
					else
						lsb_const_bits++;

					if(prev != nullptr)
						break;
				}
				continue;
			}

			SigSpec prev_B = sigmap(prev->getPort(ID::B));
			bool B_is_negative = (prev->type == ID($sub) ^ block_offset_neg);

			if(prev_B.is_fully_const()) {
				// usually B is constant and we follow A
				block_offset += (B_is_negative) ? -prev_B.as_const().as_int() : prev_B.as_const().as_int();
				prev = nullptr;
				for(auto bit : prev_A.bits()) {
					if(bit.is_wire())
						prev = sigbit_to_driver.at(bit, nullptr);
					else
						lsb_const_bits++;

					if(prev != nullptr)
						break;
				}
			} else if(prev_A.is_fully_const()) {
				// otherwise try to follow B
				block_offset += (block_offset_neg) ? -prev_A.as_const().as_int() : prev_A.as_const().as_int();
				prev = nullptr;
				for(auto bit : sigmap(prev_B).bits()) {
					if(bit.is_wire())
						prev = sigbit_to_driver.at(bit, nullptr);
					else
						lsb_const_bits++;

					if(prev != nullptr)
						break;
				}
			} else {
				// none is constant and we stop trying to simplify
				// it may still be optimizable with more effort
				break;
			}
		}

		if(verbose &&  prev) log("  %s: stopped following at %s\n", log_id(module), log_id(prev));
		if(verbose && !prev) log("  %s: stopped following, no ending cell found\n", log_id(module));

		shiftxnode.block_offset = block_offset;

		// after the block_offset we expect a multiplication (defining the block-size)
		if(prev == nullptr || prev->type != ID($mul))
			return;

		SigSpec mul_A = sigmap(prev->getPort(ID::A));
		SigSpec mul_B = sigmap(prev->getPort(ID::B));
		bool mul_B_const = false;

		if(mul_B.is_fully_const()) {
			mul_B_const = true;
			block_size = mul_B.as_const().as_int();
		} else if(mul_A.is_fully_const()) {
			block_size = mul_A.as_const().as_int();
		}

		block_size <<= lsb_const_bits;
		shiftxnode.block_size = block_size;
		shiftxnode.shift_left = block_offset_neg;
		shiftxnode.select_idx = mul_B_const ? mul_A : mul_B;

		// if one is constant, we have a shiftxnode with shift-amount = idx*block_size + block_offset
		// check the arguments
		if(block_size < min_shift)
			return;

		if(napot && ((block_size & (block_size - 1)) == 0))
			return;

		candiadates.push_back(shiftxnode);

	}


	void replace_shiftx(shiftxnode_t node) {
		SigSpec data = sigmap(node.shiftx->getPort(ID::A));
		int total_data_size = 1 << ceil_log2(node.block_offset + data.size());
		// SigSpec(const std::vector<RTLIL::SigChunk> &chunks);
		State fillbit = (node.shiftx->type == ID($shiftx)) ? State::Sx : State::S0;
		std::vector<SigChunk> chunks;

		if(node.shift_left == false) {
			// shift-right, block-mux:
			// 1. res = [ 9 8 7 6 5 4 3 2 1 0 ] >> (idx*block_size + block_offset);
			// 2. res = [ fill 9 8 7 6 5 4 3 2 1 0 ] >>  idx*block_size;
			// 3. [ fill 9 | 8 7 6  | 5 4 3 | 2 1 0 ] <-> [ B3 ] [ B2 ] [ B1 ] [ B0 ]
			// 4. res = idx >=2 ? ( idx>2 ? B3 : B2 ) : ( idx>0 ? B1: B0 );
			//addBmux (RTLIL::IdString name, const RTLIL::SigSpec &sig_a, const RTLIL::SigSpec &sig_s, const RTLIL::SigSpec &sig_y, const std::string &src = "");

			for(auto chunk: data.chunks()) // data
				chunks.push_back(chunk);
			chunks.push_back(SigChunk(fillbit, node.block_offset)); // offset
			chunks.push_back(SigChunk(fillbit, total_data_size - (node.block_offset + data.size()))); // padding to POW2

			// build index signal
			// connect to data wire
			// disconnect shiftx output

			if(verbose) log("  %s: transforming %s into $bmux: ", log_id(module), log_id(node.shiftx));
			if(verbose) log("sel-in: %s; data-in: %s; offset: %d, block-size: %d)\n",
								log_signal(node.select_idx), log_signal(data), node.block_offset, node.block_size);
			Cell *bmux = module->addBmux(NEW_ID, nullptr, nullptr, nullptr, "");
			// create bmux
			// connect bmux
			// disconnect shift/shiftx
		} else {
			// shift-left, demux:
			// this is used if something is written to only one LHS slice (and the rest is filled with whatever)
			// meaning LHS is wider than RHS and the RHS is shifted to the right position before being masked and assigned to LHS
			// 1. res =        [ 9 8 7 6 5 4 3 2 1 0 ] << (idx*block_size + block_offset);
			// 2. res = [ 9 8 7 6 5 4 3 2 1 0 block_offset ] <<  idx*block_size;
			// 3. [ 9 8 7 | 6 5 4 | 3 2 1 | 0 block_offset ]

			chunks.push_back(SigChunk(fillbit, node.block_offset)); // offset
			for(auto chunk: data.chunks()) // data
				chunks.push_back(chunk);
			chunks.push_back(SigChunk(fillbit, total_data_size - (node.block_offset + data.size()))); // padding to POW2

			if(verbose) log("  %s: transforming %s into $demux: ", log_id(module), log_id(node.shiftx));
			if(verbose) log("sel-in: %s; data-in: %s; offset: %d, block-size: %d)\n",
								log_signal(node.select_idx), log_signal(data), node.block_offset, node.block_size);
		}
	}

	void run()
	{
		log("Running $shift/$shiftx optimization on module %s:\n", log_id(module));

		index();

		for(auto cell : shiftx_cells)
			check_shiftx(cell);

		for(auto node : candiadates)
			replace_shiftx(node);

	}
};

struct OptShiftxPass : public Pass {
	OptShiftxPass() : Pass("opt_shiftx", "transform $shift/$shiftx cells, shifting by constant-sized bit-groups, to $bmux cells") { }
	void help() override
	{
		//   |---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|
		log("\n");
		log("    opt_shiftx [options] [selection]\n");
		log("\n");
		log("This pass transforms suitable $shiftx/$shiftx cells to $bmux cells.\n");
		log("$shiftx are for example used to represent Verilog part-selects, such as:\n");
		log("    out = data[15*idx:+7]; // implemented by shifting data by 15*idx.\n");
		log("Using a general shifter is often less efficient than using block-muxes ($bmux).\n");
		log("\n");
		log("    -min_shift <int>\n");
		log("        the minimum step-size being considered for this transform\n");
		log("        default: 2\n");
		log("\n");
		log("    -napot\n");
		log("        turn-on convertion of natural-power-of-two shifts to $bmux\n");
		log("\n");
		log("    -v\n");
		log("        verbose output\n");
	}
	void execute(std::vector<std::string> args, RTLIL::Design *design) override
	{
		int min_shift = 2;
		bool verbose = false;
		bool napot = false;

		log_header(design, "Executing OPT_SHIFTX pass.\n");

		size_t argidx;
		for (argidx = 1; argidx < args.size(); argidx++) {
			if (args[argidx] == "-min_shift" && argidx+1 < args.size()) {
				min_shift = atoi(args[++argidx].c_str());
				continue;
			}
			if (args[argidx] == "-napot") {
				napot = true;
				continue;
			}
			if (args[argidx] == "-v") {
				verbose = true;
				continue;
			}
			break;
		}
		extra_args(args, argidx, design);

		for (auto module : design->selected_modules())
		{
			if (!module->has_processes_warn()) {
				OptShiftxWorker worker(design, module, min_shift, napot, verbose);
				worker.run();
			}
		}

	}
} OptShiftxPass;

PRIVATE_NAMESPACE_END

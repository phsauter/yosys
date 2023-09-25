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

USING_YOSYS_NAMESPACE
PRIVATE_NAMESPACE_BEGIN

bool did_something;

#include "passes/pmgen/peepopt_pm.h"
#include "generate.h"

struct PeepoptPass : public Pass {
	PeepoptPass() : Pass("peepopt", "collection of peephole optimizers") { }
	void help() override
	{
		//   |---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|
		log("\n");
		log("    peepopt [options] [selection]\n");
		log("\n");
		log("This pass applies a collection of peephole optimizers to the current design.\n");
		log("\n");
	}
	void execute(std::vector<std::string> args, RTLIL::Design *design) override
	{
		std::string genmode;

		log_header(design, "Executing PEEPOPT pass (run peephole optimizers).\n");

		bool bmux_mode = false, wshift_mode = false;
		size_t argidx;
		for (argidx = 1; argidx < args.size(); argidx++)
		{
			if (args[argidx] == "-generate" && argidx+1 < args.size()) {
				genmode = args[++argidx];
				continue;
			}
			if (args[argidx] == "-bmux") {
				bmux_mode = true;
				continue;
			}
			if (args[argidx] == "-wshift") {
				wshift_mode = true;
				continue;
			}
			break;
		}
		extra_args(args, argidx, design);

		if (!genmode.empty())
		{
			if (genmode == "shiftmul")
				GENERATE_PATTERN(peepopt_pm, shiftmul_right);
			else if (genmode == "muldiv")
				GENERATE_PATTERN(peepopt_pm, muldiv);
			else
				log_abort();
			return;
		}

		for (auto module : design->selected_modules())
		{
			did_something = true;

			while (did_something)
			{
				did_something = false;

				peepopt_pm pm(module);

				pm.setup(module->selected_cells());

				pm.ud_shiftmul_right.bmux_mode = bmux_mode;
				pm.ud_shiftmul_right.wshift_mode = wshift_mode;

				pm.run_shiftmul_right();
				pm.run_shiftmul_left();
				pm.run_muldiv();
			}
		}
	}
} PeepoptPass;

PRIVATE_NAMESPACE_END

import asyncio
from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, ResultMessage

WORKING_DIR = "/Users/va/Documents/GitHub/Vasilis local/Claude Agent 2nd try"
MAX_LOOPS = 2

async def run_agent(prompt, label, model="claude-opus-5", effort="xhigh"):
    """Run one agent, print its output, and return the full response as text."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}\n")
    
    response = []
    async for message in query(
        prompt=prompt,
        options=ClaudeAgentOptions(
            allowed_tools=["Read", "Write", "Bash"],
            permission_mode="acceptEdits",
            cwd=WORKING_DIR,
            model=model,
            effort=effort,              # xhigh suits coding/agentic work on Opus 5
            max_turns=40,               # ceiling so a confused agent can't loop forever
            max_budget_usd=5.00,        # hard dollar cap per agent run
            fallback_model="claude-opus-4-8",
            # A single SDK message must hold the whole tool result. Reading
            # thz_results.png (705 KB -> 940 KB base64) blew the 1 MB default
            # and killed the run after the coder had finished its work.
            max_buffer_size=32 * 1024 * 1024,
        ),
    ):
        if isinstance(message, AssistantMessage):
            for block in message.content:
                if hasattr(block, "text"):
                    print(block.text, end="", flush=True)
                    response.append(block.text)
        elif isinstance(message, ResultMessage):
            print(f"\n[{label} done: {message.subtype}]")
    
    return "".join(response)


async def main():
    reviewer_feedback = ""  # starts empty, fills after first review

    for loop in range(1, MAX_LOOPS + 1):
        print(f"\n\n{'#'*60}")
        print(f"  LOOP {loop} of {MAX_LOOPS}")
        print(f"{'#'*60}")

        # ── AGENT A: Coder ──────────────────────────────────────
        if loop == 1:
            coder_prompt = """
You are a THz-TDS spectroscopy expert and Python developer.

This folder contains two data files from a Menlo Systems THz-TDS setup:
- A sample measurement (Si wafer ~385 µm on a quartz window ~1 mm)
- A reference measurement (empty beam path)

Explore the folder, find the files, and write a script called 'thz_analysis.py' that:
1. Reads both data files
2. FFTs both waveforms (zero-pad to 2^14 points, apply a Hann window)
3. Computes complex transmission T(ω) = E_sample / E_ref
4. Uses a Fabry-Perot transfer matrix model for the Si/quartz bilayer
5. Extracts n(ω) and k(ω) for each layer separately
6. Estimates layer thickness d from echo timing in the time domain
7. Saves a clear plot as 'thz_results.png' with labelled axes and units
8. Prints a summary table of n, k, d at key frequencies

Then run the script and report the results.
Write well-commented, physically correct code.
"""
        else:
            coder_prompt = f"""
You are a THz-TDS spectroscopy expert and Python developer.

You previously wrote 'thz_analysis.py'. The reviewer has provided this feedback:

{reviewer_feedback}

Update 'thz_analysis.py' to address every point raised.
Then run the updated script and report the new results.
"""

        coder_output = await run_agent(coder_prompt, f"AGENT A — Coder (loop {loop})")

        # ── AGENT B: Reviewer ────────────────────────────────────
        reviewer_prompt = f"""
You are an expert in THz time-domain spectroscopy and scientific Python.

The coder has just written or updated 'thz_analysis.py'. Here is their report:

{coder_output}

Your job:
1. Read 'thz_analysis.py' carefully
2. Check the physics: is the Fabry-Perot model correct? Are n, k, d physically reasonable for Si and quartz?
3. Check numerical stability and edge cases (division by zero near DC, noise floor handling, etc.)
4. Check the plot: axes labelled? units shown? frequency range sensible?
5. Check code clarity and comments

If you find meaningful improvements, list them clearly and specifically.

If the code is correct and complete, respond with exactly:
NO FURTHER IMPROVEMENTS NEEDED
...followed by a brief summary of the final results.
"""

        reviewer_output = await run_agent(reviewer_prompt, f"AGENT B — Reviewer (loop {loop})")

        # ── Check stopping condition ─────────────────────────────
        # Anchor on the FIRST line only: the sentinel phrase also appears in the
        # reviewer prompt, which is echoed through coder_output, so a substring
        # search anywhere in the text can end the loop by accident.
        first_line = next(
            (ln.strip() for ln in reviewer_output.splitlines() if ln.strip()), ""
        )
        if first_line.upper().startswith("NO FURTHER IMPROVEMENTS NEEDED"):
            print("\n✅  Reviewer satisfied — pipeline complete.")
            break

        reviewer_feedback = reviewer_output  # pass to next coder loop

        if loop == MAX_LOOPS:
            print(f"\n⏹  Reached maximum of {MAX_LOOPS} loops.")

asyncio.run(main())
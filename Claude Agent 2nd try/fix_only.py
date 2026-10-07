"""Coder-only run: Agent A applying Agent B's review findings to thz_analysis.py.

This is the loop-2 coder step from 'agent 2.py', run standalone. The original
pipeline crashed in loop 1 before the reviewer ever ran, so the review was done
separately by 'review_only.py'; its output is in 'review.log' and gets fed in
here as reviewer_feedback.

Pre-fix copies of thz_analysis.py, thz_results.png and review.log are in
pre_fix_backup/ .
"""

import asyncio
from pathlib import Path

from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, ResultMessage

WORKING_DIR = "/Users/va/Documents/GitHub/Vasilis local/Claude Agent 2nd try"

# review.log opens with the reviewer's live tool-use chatter; the findings proper
# begin at this line. Slice from there so the coder gets the report, not the noise.
_raw = Path(WORKING_DIR, "review.log").read_text(encoding="utf-8")
_marker = "I verified everything by running the code"
_idx = _raw.find(_marker)
if _idx == -1:
    raise SystemExit("could not locate the findings section in review.log")
reviewer_feedback = _raw[_idx:].strip()


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
            effort=effort,
            max_turns=60,  # more than the reviewer: this one edits, runs, re-checks
            max_budget_usd=8.00,
            fallback_model="claude-opus-4-8",
            max_buffer_size=32 * 1024 * 1024,  # the PNG read that killed the first run
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


coder_prompt = f"""
You are a THz-TDS spectroscopy expert and Python developer.

You previously wrote 'thz_analysis.py'. The reviewer has provided this feedback:

{reviewer_feedback}

Update 'thz_analysis.py' to address every point raised.
Then run the updated script and report the new results.

Additional requirements for this pass:

- Finding 1 (group vs phase delay) and finding 2 (the Si upper-bound k leaking
  into the quartz inversion) are real physics bugs. Fix those first and state
  explicitly how each headline number moves as a result.

- N_SI_REF = 3.4175 on line 170 is an ASSUMED external anchor, not a measurement,
  and every quartz number inherits it. Do not present n_Si as a result of the
  analysis. Label it clearly as an input wherever it is printed or plotted, and
  make the anchor's contribution to the quartz uncertainty explicit.

- The previous run reported quartz alpha = 0.31 / 2.6 / 8.2 cm^-1, but the code
  on disk prints 0.157 / 2.11 / 7.82. Make sure the numbers you report are the
  numbers the committed script actually produces. Run it and quote its output.

- Quote results with uncertainties at a defensible precision. Do not clip
  negative absorption silently; an unresolved bin should read as unresolved,
  not as zero.

- IMPORTANT -- the reviewer is wrong about nominal thickness, and you must not
  inherit that error. The layer thicknesses are NOT known. d is one of the
  parameters this analysis EXTRACTS. The quoted "~385 um" and "~1 mm" are loose
  manufacturer figures with their own uncertainty, not reference values.
  Specifically reject these two pieces of the review:
    * Its argument that the S_main phase estimator is better because the result
      "moves toward literature fused silica and the 1000 um nominal". Choose the
      estimator on physics and numerics alone -- which one is free of the
      Fabry-Perot fringe phase, which is better conditioned, which has smaller
      frequency-resolved scatter. Proximity to a nominal number is not evidence.
    * Its claim that "d_Si = 385.4 um vs 385 nominal (+0.10 %)" is an impressive
      independent check. It is not independent: d_Si = (n*d)_Si / N_SI_REF, so it
      tests the assumed anchor against a loose nominal figure, and neither side
      is ground truth. Remove or correct that framing wherever it appears.
  In the output, state plainly which quantities are actually measured (the
  optical thickness n*d), which are assumed inputs (the index anchor), and which
  are derived from splitting the two (n and d separately, strongly
  anti-correlated). Do not present agreement with a nominal value as validation.

If you disagree with any finding, say so and justify it rather than changing the
code -- a wrong "fix" is worse than a documented disagreement.
"""


async def main():
    await run_agent(coder_prompt, "AGENT A — Coder (applying review)")


asyncio.run(main())

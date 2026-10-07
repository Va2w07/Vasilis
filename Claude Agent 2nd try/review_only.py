"""Review-only run: Agent B alone, against the thz_analysis.py already on disk.

The full pipeline in 'agent 2.py' crashed during Loop 1 / Agent A when the SDK
tried to read back the 705 KB thz_results.png (940 KB base64 > the 1 MB default
max_buffer_size). The coder had already finished its work, so the script and
figure are good -- what never happened was the review. This runs just that.
"""

import asyncio
from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, ResultMessage

WORKING_DIR = "/Users/va/Documents/GitHub/Vasilis local/Claude Agent 2nd try"


async def run_agent(prompt, label, model="claude-opus-5", effort="xhigh"):
    """Run one agent, print its output, and return the full response as text."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}\n")

    response = []
    async for message in query(
        prompt=prompt,
        options=ClaudeAgentOptions(
            allowed_tools=["Read", "Bash"],  # reviewer reads and may re-run; it does not edit
            permission_mode="acceptEdits",
            cwd=WORKING_DIR,
            model=model,
            effort=effort,
            max_turns=40,
            max_budget_usd=5.00,
            fallback_model="claude-opus-4-8",
            max_buffer_size=32 * 1024 * 1024,  # the fix: PNG reads blew the 1 MB default
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


# Stands in for {coder_output} in the original prompt -- the coder's own report was
# lost with the crash, so point the reviewer at the artefacts instead.
CODER_REPORT = """
(The coder's live report was lost when the run crashed. Its work is on disk:
'thz_analysis.py' and 'thz_results.png'. From its partial log, it reported:

- 2425 points, 80 ps window, 33 fs steps; sample pulse delayed 6.30 ps
- Echoes at +8.78 ps (Si round trip), +13.00 ps (quartz), +21.8 ps (cross term)
- Found system ringing in the REFERENCE at +7.23/+7.89/+8.94 ps landing on the
  Si echo gate; worked around it by gating the deconvolved impulse response
  rather than the raw waveform
- Final: n_Si = 3.4175, n_quartz = 1.9762, d_Si = 385.4 um, d_quartz = 988.3 um,
  quartz alpha = 0.31 / 2.6 / 8.2 cm^-1 at 0.5 / 1.0 / 2.0 THz,
  Si alpha <= 1.2 cm^-1 upper bound, transfer-matrix residual 0.0098

Verify these claims against the code yourself rather than taking them on trust.)
"""

reviewer_prompt = f"""
You are an expert in THz time-domain spectroscopy and scientific Python.

The coder has just written or updated 'thz_analysis.py'. Here is their report:

{CODER_REPORT}

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


async def main():
    await run_agent(reviewer_prompt, "AGENT B — Reviewer (review-only run)")


asyncio.run(main())

import asyncio
from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, ResultMessage

async def main():
    async for message in query(
        prompt="""
These are THz time-domain spectroscopy files from a Menlo Systems setup.
The sample is a Si wafer (~385 µm) on a quartz window (~1 mm).
The reference file is the empty beam path.

Write a Python script called 'thz_analysis.py' in this folder that:
1. Reads both data files
2. FFTs both waveforms with zero-padding
3. Computes the complex transmission T(ω) = E_sample / E_ref
4. Uses a Fabry-Perot transfer matrix model for the Si/quartz stack
5. Extracts n(ω) and k(ω) for each layer separately
6. Estimates d for each layer from the echo timing
7. Plots the results and saves the figures as PNG files

Use numpy, scipy, and matplotlib. Add comments explaining each step.
""",
        options=ClaudeAgentOptions(
            allowed_tools=["Read", "Bash"],
            permission_mode="acceptEdits",
            cwd="/Users/va/Documents/GitHub/Vasilis local/Claude Agent 1st try",
        ),
    ):
        if isinstance(message, AssistantMessage):
            for block in message.content:
                if hasattr(block, "text"):
                    print(block.text)
        elif isinstance(message, ResultMessage):
            print(f"Done: {message.subtype}")

asyncio.run(main())
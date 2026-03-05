import asyncio
from backend.src.api.conditions import _run_conditions, ConditionGridRequest
from backend.src.utils.config import load_config

async def main():
    config = load_config()
    request = ConditionGridRequest(
        max_dist=[10, 15],
        compactness=[0.06, 0.1],
        dilation=[5],
        roi_name="text"
    )
    await _run_conditions(config, request)

asyncio.run(main())

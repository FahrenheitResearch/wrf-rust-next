# GPU Experimental Path

Use this only when the user explicitly wants GPU WRF or is already stuck on GPU-specific crashes.

## Default stance

Do not put beginners here first.

Require:

- a working CPU baseline
- clear user willingness to debug compiler and runtime issues
- realistic expectations

## What local evidence shows

- A `3 km` GPU run crashed with `CUDA_ERROR_ILLEGAL_ADDRESS` in `module_bl_ysu.F`.
- A `1 km`, `400 x 400 x 80` case used about `31.2 GiB` on a `32 GiB` GPU.
- There were GPU-only device-sync issues around output handling.

## How to message it

"The GPU branch is experimental. Get the CPU path stable first, then try the GPU build."

## Fast troubleshooting rule

If the user reports:

- `CUDA_ERROR_ILLEGAL_ADDRESS`
- `module_bl_ysu.F`
- device data missing
- weird OpenACC runtime faults

tell them to:

1. confirm the CPU run works
2. test a smaller case
3. try a host fallback or simplified physics branch
4. stop debugging GPU if the real goal is just getting a forecast out

## When to advise against it

- Windows user with no working WSL CPU baseline
- 16 to 32 GB host RAM and tiny patience
- user only wants a stable first run


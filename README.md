# TetherModels

[![Build Status](https://github.com/AndreaBertozzi/TetherModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/AndreaBertozzi/TetherModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Goals
1. provide a quasi-static tether model for realistic testing of [WinchControllers.jl](https://github.com/OpenSourceAWE/WinchControllers.jl)
2. provide a more generic tether model where both end points can freely be defined for simulating e.g. a kite, pulling a boat (the lower end of the tether is moving)

**Non-goal:** This package shall NOT use modelling toolkit (MTK), at least not as obligatory dependency. Reason: MTK is very large and loading it takes a lot of time. For now, this package should be small and fast to load, so that it can well be used as test dependency for example for `WinchControllers.jl`. 

In the beginning, the focus should be on implementing the first goal. Working towards the second goal can come later.
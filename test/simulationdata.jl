using DynamicGrids, Test

init = [1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1]

ruleset = Ruleset()
starttime = 1
radius = 3

simdata = DynamicGrids.SimData(init, ruleset, starttime)

# TODO some tests
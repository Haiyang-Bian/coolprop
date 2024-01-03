using Oxygen, HTTP

include("Caculators.jl")


@post "/table/sat" function (req::HTTP.Request)
    request = json(req, TableSatArgs)
    return caculate(request)
end

@post "/table/supercooling" function (req::HTTP.Request)
    request = json(req, TableSuperCoolingArgs)
    return caculate(request)
end

@post "/table/overheated" function (req::HTTP.Request)
    request = json(req, TableOverHeatedArgs)
    return caculate(request)
end

@post "/plot/hsdiagram" function (req::HTTP.Request)
    request = json(req, HSDiagram)
    return diagram(request)
end

serve(async=true)
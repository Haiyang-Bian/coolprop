"""
网页小项目:基于CoolProp的热力学计算后端,v1.0
预计实现的基本功能:
1.纯净物的物性查询
2.混合物的物性查询
3.物性查询的基本方面:
	1.基本的热力学函数(比量):比焓kJ/kg,比熵kJ/(kg*K),密度(kg/m^3)/比体积(m^3/kg),比吉布斯自由能(kJ/kg)
	2.饱和态,过冷态,过热态.(不考虑固相,也即主要为气-液二相)
	3.临界,及临界以上
4.查询命令及参数
	targ -h/s/g/v/r (比焓/比熵/比吉布斯自由能/比体积/密度)以质量为比
	mode --table
	thermalstate --sat/overheated/supercool/overcritical
	T 温度(绝对温标/K)
	P 压强(可指定单位,默认为Pa)
	Q 干度(只有饱和态才可使用)
	example: coolprop water targ -h T 300K P 50.2MPa 
	/table
	{
		"media":"water",
		"thermalstate":"sat",
		"orderedby":"T",
		"begin":"200K",
		"end":"500K",
		"step":"5K",
		"targ":["H","S"]
	}
	/one
	{
		"media":"water",
		"thermalstate":"sat",
		"
		"targ":["H","S"]
	}
"""

using CoolProp, Unitful, Plots

"""
饱和态物性表查询信息体.
以压强为序时,步长单位默认为MPa
"""
struct TableSatArgs
    media::String
    targ::Vector{String}
    orderedby::String
    start::String
    ends::String
    step::String
end
"""
饱和态物性查询计算函数.
"""
function caculate(input::TableSatArgs)
    b = uparse(input.start)
    e = uparse(input.ends)
    step = uparse(input.step)
    e < b && error("终止位置要大于起始位置!")
    if input.orderedby == "T"
        tcrit = PropsSI("Tcrit", input.media) * u"K"
        t0 = PropsSI("Ttriple", input.media) * u"K"
        tspan = max(t0, b):input.step:min(e, tcrit)
        pspan = PropsSI.("P", "T", tspan, "Q", 0.0, input.media)
        ans = Dict("T" => tspan, "P" => pspan)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "T", tspan, "Q", 0.0, input.media))
            get!(ans, "$tag-Gas", PropsSI.(tag, "T", tspan, "Q", 1.0, input.media))
        end
        return ans
    elseif input.orderedby == "P"
        pcrit = PropsSI("Pcrit", input.media) * u"Pa"
        p0 = PropsSI("P_min", input.media) * u"Pa"
        pspan = max(b, p0):step:min(pcrit, e)
        tspan = PropsSI.("T", "P", pspan, "Q", 0.0, input.media)
        ans = Dict("P" => pspan, "T" => tspan)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "P", pspan, "Q", 0.0, input.media))
            get!(ans, "$tag-Gas", PropsSI.(tag, "P", pspan, "Q", 1.0, input.media))
        end
        return ans
    else
        error("错误!排序只能以温度或压强!")
    end
end
"""
过冷液体物性计算信息体
以压强为序时,单位为MPa.
"""
struct TableSuperCoolingArgs
    media::String
    targ::Vector{String}
    orderedby::String
    start::String
    ends::String
    step::String
    given::String
end
"""
过冷液体物性计算函数
"""
function caculate(input::TableSuperCoolingArgs)
    given = uparse(input.given)
    step = uparse(input.step)
    start = uparse(input.start)
    ends = uparse(input.ends)
    if input.orderedby == "T"
        t = PropsSI("T", "P", given, "Q", 0.0, input.media)
        t0 = PropsSI("T_min", input.media) * u"K"
        tspan = max(t0, start):step:min(t, ends)
        ans = Dict("P" => given)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "T", tspan, "P", given, input.media))
        end
        return ans
    elseif input.orderedby == "P"
        p = PropsSI("P", "T", given, "Q", 0.0, input.media)
        p0 = PropsSI("P_min", input.media) * u"Pa"
        pspan = max(p0, start):step:min(p, ends)
        ans = Dict("T" => given)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "P", pspan, "T", given, input.media))
        end
        return ans
    else
        error("出错了!")
    end
end
"""
过热蒸汽物性计算信息体(其结构与过冷相似)
以压强为序时,单位为MPa.
"""
struct TableOverHeatedArgs
    media::String
    targ::Vector{String}
    orderedby::String
    start::String
    ends::String
    step::String
    given::String
end
"""
过热蒸汽物性计算函数
"""
function caculate(input::TableOverHeatedArgs)
    given = uparse(input.given)
    step = uparse(input.step)
    start = uparse(input.start)
    ends = uparse(input.ends)
    if input.orderedby == "T"
        t = PropsSI("Tcrit", input.media) * u"K"
        t0 = PropsSI("T", "P", given, "Q", 1.0, input.media)
        tspan = max(t0, start):step:min(t, ends)
        ans = Dict("P" => given)
        for tag in input.targ
            get!(ans, "$tag-Gas", PropsSI.(tag, "T", tspan, "P", given, input.media))
        end
        return ans
    elseif input.orderedby == "P"
        p = PropsSI("Pcrit", input.media) * u"Pa"
        p0 = PropsSI("P", "T", given, "Q", 1.0, input.media)
        pspan = max(p0, start):step:min(p, ends)
        ans = Dict("T" => given)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "P", pspan, "T", given, input.media))
        end
        return ans
    else
        error("出错了!")
    end
end
"""
超临界物性计算信息体(其结构与过冷相似)
以压强为序时,单位为MPa.
"""
struct TableOverCriticalArgs
    media::String
    targ::Vector{String}
    orderedby::String
    start::String
    ends::String
    step::String
    given::String
end
"""
超临界物性计算函数
"""
function caculate(input::TableOverCriticalArgs)
    given = uparse(input.given)
    step = uparse(input.step)
    start = uparse(input.start)
    ends = uparse(input.ends)
    if input.orderedby == "T"
        t0 = PropsSI("Tcrit", input.media) * u"K"
        t = PropsSI("T_max", input.media) * u"K"
        tspan = max(t0, start):step:min(t, ends)
        ans = Dict("P" => given)
        for tag in input.targ
            get!(ans, "$tag-Gas", PropsSI.(tag, "T", tspan, "P", given, input.media))
        end
        return ans
    elseif input.orderedby == "P"
        p0 = PropsSI("Pcrit", input.media) * u"Pa"
        p = PropsSI("P_max", input.media) * u"Pa"
        pspan = max(p0, start):step:min(p, ends)
        ans = Dict("T" => given)
        for tag in input.targ
            get!(ans, "$tag-Liquid", PropsSI.(tag, "P", pspan, "T", given, input.media))
        end
        return ans
    else
        error("出错了!")
    end
end
"""
起始,步长,终止信息体
"""
struct StartStepStop
    start::String
    step::String
    stop::String
end
"""
重载迭代器
"""
Base.StepRange(sss::StartStepStop) = StepRange(uparse(sss.start), uparse(sss.step), uparse(sss.stop))
global points = 50
"""
绘制焓熵图的信息体
"""
struct HSDiagram
    media::String
    tspan::StartStepStop
    pspan::StartStepStop
    vspan::StartStepStop
    qspan::StartStepStop
end
"""
绘制焓熵图的函数
"""
function diagram(input::HSDiagram)
    ans = Dict()
    points_num = points
    e = 0.1
    if !isnothing(input.tspan)
        # 绘制等温线
        p0 = (PropsSI("P_min", input.media) + e) * u"Pa"
        pm = (PropsSI("P_max", input.media) - e) * u"Pa"
        pspan = range(p0, pm, points_num)
        t_dict = Dict()
        for line in StepRange(input.tspan)
            s = PropsSI.("S", "P", pspan, "T", line, input.media)
            h = PropsSI.("H", "P", pspan, "T", line, input.media)
            line_dict = Dict("s" => s, "h" => h)
            get!(t_dict, line, line_dict)
        end
        get!(ans, "T", t_dict)
    end
    if !isnothing(input.pspan)
        # 绘制等压线
        t0 = (PropsSI("T_min", input.media) + e) * u"K"
        tm = (PropsSI("T_max", input.media) - e) * u"K"
        tspan = range(t0, tm, points_num)
        p_dict = Dict()
        for line in StepRange(input.pspan)
            s = PropsSI.("S", "T", tspan, "P", line, input.media)
            h = PropsSI.("H", "T", tspan, "P", line, input.media)
            line_dict = Dict("s" => s, "h" => h)
            get!(p_dict, line, line_dict)
        end
        get!(ans, "P", p_dict)
    end
    if !isnothing(input.vspan)
        # 绘制等体积线
        t0 = (PropsSI("T_min", input.media) + e) * u"K"
        tm = (PropsSI("T_max", input.media) - e) * u"K"
        tspan = range(t0, tm, points_num)
        v_dict = Dict()
        for line in 1 ./ StepRange(input.vspan)
            s = PropsSI.("S", "T", tspan, "D", line, input.media)
            h = PropsSI.("H", "T", tspan, "D", line, input.media)
            line_dict = Dict("s" => s, "h" => h)
            get!(v_dict, line, line_dict)
        end
        get!(ans, "V", v_dict)
    end
    if !isnothing(input.qspan)
        # 绘制等干度线
        t0 = (PropsSI("T_min", input.media) + e) * u"K"
        tm = (PropsSI("Tcrit", input.media) - e) * u"K"
        tspan = range(t0, tm, points_num)
        q_dict = Dict()
        for line in StepRange(input.qspan)
            s = PropsSI.("S", "T", tspan, "Q", line, input.media)
            h = PropsSI.("H", "T", tspan, "Q", line, input.media)
            line_dict = Dict("s" => s, "h" => h)
            get!(q_dict, line, line_dict)
        end
        get!(ans, "Q", q_dict)
    end
    return ans
end
"""
绘制对数压焓图的信息体
"""
struct LogPHDiagram
    media::String
    tspan::StartStepStop
    sspan::StartStepStop
    vspan::StartStepStop
    qspan::StartStepStop
end
"""
绘制对数压焓图图的函数
"""
function diagram(input::LogPHDiagram)
    points_num = points
    e = 0.1
    p0 = (PropsSI("P_min", input.media) + e) * u"Pa"
    pm = (PropsSI("P_max", input.media) - e) * u"Pa"
    pspan = range(p0, pm, points_num)
    logp = log10.(pspan)
    ans = Dict("LogP" => logp)
    if !isnothing(input.tspan)
        # 绘制等温线
        t_dict = Dict()
        for line in StepRange(input.tspan)
            h = PropsSI.("H", "P", pspan, "T", line, input.media)
            get!(t_dict, line, h)
        end
        get!(ans, "T", t_dict)
    end
    if !isnothing(input.sspan)
        # 绘制等熵线
        s_dict = Dict()
        for line in StepRange(input.sspan)
            h = PropsSI.("H", "P", pspan, "S", line, input.media)
            get!(s_dict, line, h)
        end
        get!(ans, "S", s_dict)
    end
    if !isnothing(input.vspan)
        # 绘制等体积线
        v_dict = Dict()
        for line in 1 ./ StepRange(input.vspan)
            h = PropsSI.("H", "P", pspan, "D", line, input.media)
            get!(v_dict, line, h)
        end
        get!(ans, "V", v_dict)
    end
    if !isnothing(input.qspan)
        # 绘制等干度线
        q_dict = Dict()
        for line in StepRange(input.qspan)
            h = PropsSI.("H", "P", pspan, "Q", line, input.media)
            get!(q_dict, line, h)
        end
        get!(ans, "Q", q_dict)
    end
    return ans
end

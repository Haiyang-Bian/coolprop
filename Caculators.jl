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

using CoolProp, Unitful


struct Args
    media
    targ::String
    mode
    thermalstate
    t
    p
    q
end
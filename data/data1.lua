numlab=require"numlab"

function printf(...)
	io.write(string.format(...))
end

local w=logspace(10,1000,100)
local w0=100
local Ta=abp2(w,w0,0.25)
local T=bp2(w,w0,0.25)

for i=1,#w do
	printf("%g %g\n",w[i],(abs(T))[i])
end

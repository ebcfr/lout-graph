numlab=require"numlab"

function printf(...)
	io.write(string.format(...))
end

E=1;f0=1e3;m=0.1;T0=1
w0=2*pi*f0
t=linspace(-1e-3,4e-3,200)
u=E*ustep(t)
y=lsim({T0},{1,2*m/w0,1/w0^2},u,t)
y1=T0*E*(1+1/math.sqrt(1-m^2)*exp(-m*w0*t))*ustep(t)
y2=T0*E*(1-1/math.sqrt(1-m^2)*exp(-m*w0*t))*ustep(t)

dt=t[2]-t[1]
t1=0*1e-3;t2=3e-3
i1=math.floor((t1+1e-3)/dt)
i2=math.floor((t2+1e-3)/dt)

f=io.open("data2","w")
for i=1,#t do
  f:write(t[i]*1000 .. " " .. y[1][i] .. "\n")
end
f:close()

f=io.open("data3","w")
for i=1,#t do
  f:write(t[i]*1000 .. " " .. y1[i] .. "\n")
end
f:close()

f=io.open("data4","w")
for i=1,#t do
  f:write(t[i]*1000 .. " " .. y2[i] .. "\n")
end
f:close()

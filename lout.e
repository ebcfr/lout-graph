j=1i;

function residue(num,den,toler=0.001)
	if length(num)==0 || length(den)==0;error("a or b is not a polynom");endif;
	a=den;b=num;
	la = length(a);lb=length(b);
	if la==1;return {[],[],b/a,[]};endif;
	b=b/a[la];a=a/a[la];
	
	.. find the poles
	p=polysolve(a);lp=length(p);
	
	.. are the poles real ?
	index = nonzeros(abs(im(p)/re(p)) < toler);
	if (length(index) != 0)
		p[index] = re(p[index]);
	endif
	
	.. Find the direct term if there is one.

	if (lb >= la)	.. Also returns the reduced numerator.
		{k, b} = polydiv(b, a);
		b=polytrunc(b);
		lb = length (b);
	else
		k = [];
	endif

	if (lp == 1)
		r = polyval (b, p);
		e = 1;
		return {r,p,k,e};
	endif

	.. We need to determine the number and multiplicity of the roots.
	..
	.. D  is the number of distinct roots.
	.. M  is a vector of length D containing the multiplicity of each root.
	.. pr is a vector of length D containing only the distinct roots.
	.. e  is a vector of length lp which indicates the power in the partial
	..    fraction expansion of each term in p.

	.. Set initial values.  We'll remove elements from pr as we find
	.. multiplicities.  We'll shorten M afterwards.

	e = ones (1,lp);
	M = zeros (1,lp);
	pr = p;
	D = 1;
	M[1] = 1;

	oldpindex = 1;
	Mindex = 1;
	prindex = 2;

	loop 2 to lp;
		if (abs (p[#] - p[oldpindex]) < toler)
			.. We've found a multiple pole.
			M[Mindex] = M[Mindex] + 1;
			e[#] = e[#-1] + 1;
			.. Remove the pole from pr.
			pr[prindex:length(pr)-1] = pr[prindex+1:length(pr)];
			pr[length(pr)]=0;
		else
			.. It's a different pole.
			D=D+1;
			Mindex=Mindex+1;
			M[Mindex] = 1;
			oldpindex = #;
			prindex=prindex+1;
		endif
	end

	.. Shorten M to it's proper length

	M = M[1:D];
	
	MM = max(M);

	.. Left hand side polynomial

	lhs = zeros(MM, lb);
	rhs = zeros(MM, lp*lp);
	lhs[1,:] = b;
	rhi = 1;
	mpi = 1;
	loop 1 to D
		for ind = 1 to M[#]
			if (mpi > 1 && (mpi+ind) <= lp)
				cp = p[1:mpi-1]|p[mpi+ind:lp];
			elseif (mpi == 1)
				cp = p[mpi+ind:lp];
			else
				cp = p[1:mpi-1];
			endif
			rhs[1, rhi:rhi+lp-ind]=polycons(cp);
			rhi = rhi + lp;
		end
		mpi = mpi + M[#];
	end
	if (MM > 1)
		for index = 2 to MM
			lhs[index,1:lb-index+1]=polydif(lhs[index-1,1:lb-index+2]);
			ind = 1;
			for rhi = 1 to lp
				cp = rhs[index-1, ind:ind+lp-index+1];
				rhs[index, ind:ind+lp-index] = polydif(cp);
				ind = ind + lp;
			end
		end
	endif

	.. Now lhs contains the numerator polynomial and as many derivatives as
	.. are required.  rhs is a matrix of polynomials, the first row
	.. contains the corresponding polynomial for each residue and
	.. successive rows are derivatives.

	.. Now we need to evaluate the first row of lhs and rhs at each
	.. distinct pole value.  If there are multiple poles we will also need
	.. to evaluate the derivatives at the pole value also.

	B = zeros (lp, 1);
	A = zeros (lp, lp);

	row = 1;
	loop 1 to D
		for mi = 1 to M[#]
			B[row] = polyval (lhs[mi, :], pr[#]);
			ci = 1;
			for col = 1 to lp
				cp = rhs[mi, ci:ci+lp-1];
				A[row, col] = polyval(cp,pr[#]);
				ci = ci + lp;
			end
			row=row+1;
		end
	end

	.. Solve for the residues.

	r = A \ B;
	
	return {r',p,k,e};
endfunction

function partial(num,den,x="s")
	{r,p,k,e}=residue(num,den);
	out="";
	n=length(k);
	if n
		loop 1 to n
			if (im(k[#])~=0);
				out=out|printf("%g",re(k[#]));
			else
				out=out|" + "|printf("%g",re(k[#]));
				if im(k[#])>0
					out=out|printf("+%gj",im(k[#]));
				else
					out=out|printf("%gj",im(r[#]));
				endif
			endif
			if #>2
				out=out|" "|x|printf(" sup %g",#-1);
			elseif #==2
				out=out|" "|x;
			endif
			if #<n
				out=out|" + ";
			endif
		end
	endif;
	if length(r)
		if length(k);out=out|" + ";endif;
		loop 1 to length(r)
			signrep=(re(p[#])<0);
			if (im(r[#])~=0)
				out=out|printf("{%g}",re(r[#]))|" over {"|x;
			elseif re(r[#])~=0
				out=out|printf("{%gj}",im(r[#]))|" over {"|x;
			else
				out=out|printf("{%g",re(r[#]));
				if im(r[#])>0
					out=out|printf("+%gj",im(r[#]))|"} over {"|x;
				else
					out=out|printf("%gj",im(r[#]))|"} over {"|x;
				endif
			endif;
			if (im(p[#])~=0)
				if signrep
					out=out|printf("+ %g",-re(p[#]))|"}";
				else
					out=out|printf("- %g",re(p[#]))|"}";
				endif
			elseif re(p[#])~=0
				if im(p[#])<0
					out=out|printf("+ %g",-im(p[#]))|"}";
				else
					out=out|printf("- %g",im(p[#]))|"}";
				endif
			else
				if signrep
					out=out|printf("+ %g",-re(p[#]));
					if im(p[#])>0
						out=out|printf("- %gj",im(p[#]))|"}";
					else
						out=out|printf("+ %gj",-im(p[#]))|"}";
					endif
				else
					out=out|printf("- %g",re(p[#]));
					if im(p[#])>0
						out=out|printf("+ %gj",im(p[#]))|"}";
					else
						out=out|printf("- %gj",-im(p[#]))|"}";
					endif
				endif
			endif;
			if (e[#]>1)
				out=out|" sup "|printf("%g",e[#]);
			endif;
			if #<length(r);out=out|" + ";endif;
		end
	endif
	return out;
endfunction

function polydisp(a, x="s")
	out=""
	for i=length(a) to 1 step -1
		if iscomplex(a[i])
			if !(im(a[i])~=0) && !(re(a[i])~=0)
				if i!=length(a);" + ";endif;
				if i>2
					"(" a[i] ")`"|x|" sup "|printf("%g",i-1)
				elseif i==2
					"(" a[i] ")`"|x
				else
					if re(a[i])<0
						a[i]
					else
						" + " a[i]
					endif
				endif
			elseif im(p[i])~=0
				if i!=length(a) && re(a[i])>=0
					" + "
				endif;
				if i>2
					re(a[i]) x|" sup "|printf("%g",i-1)
				elseif i==2
					re(a[i]) x
				else
					re(a[i])
				endif
			else
				if i!=length(a) && im(a[i])>=0
					" + "
				endif;
				if i>2
					im(a[i]) x|" sup "|printf("%g",i-1)
				elseif i==2
					im(a[i]) x
				else
					im(a[i])
				endif
			endif
		else
			if i!=length(a) && a[i]>=0
				" + "
			endif;
			if i>2
				a[i] x|" sup "|printf("%g",i-1)
			elseif i==2
				a[i] x
			else
				a[i]
			endif
		endif
	end
	return out;
endfunction

function getcolor(c)
	if c==1 return "black";
	elseif c==2 return "red";
	elseif c==3 return "magenta";
	elseif c==4 return "blue";
	elseif c==5 return "cyan";
	elseif c==6 return "darkcyan";
	elseif c==7 return "darkred";
	elseif c==8 return "darkgreen";
	elseif c==9 return "darkblue";
	elseif c==10 return "darkcyan";
	elseif c==11 return "darkmagenta";
	elseif c==12 return "darkyellow";
	elseif c==13 return "lightred";
	elseif c==14 return "lightgreen";
	elseif c==15 return "lightblue";
	endif
	return "black";
endfunction

function graph()
	"@Graph "
	if isvar(width) "width{" width "} " endif
	if isvar(height) "height{" height "} " endif
	if isvar(grid) "grid{ " grid " } "  endif
	if isvar(style)
		"style{" style "} "
		if style=="frame"
			if isvar(xextra) "xextra{" xextra "} " else " xextra{0c} " endif
			if isvar(yextra) "yextra{" yextra "} " else " yextra{0c} " endif
		elseif style=="axes"
			if isvar(xorigin) "xorigin{" xorigin "} " elseif isvar(xlog) " xorigin{1} " else " xorigin{0} " endif
			if isvar(yorigin) "yorigin{" yorigin "} " elseif isvar(ylog) " yorigin{1} " else " yorigin{0} " endif
		endif
	else
		"style{frame} xextra{0c} yextra{0c} "
	endif
	if isvar(xlog) "xlog{10} " endif
	if isvar(ylog) "ylog{10} " endif
	if isvar(xticks) "xticks{" xticks "} " endif
	if isvar(yticks) "yticks{" yticks "} " endif
	if isvar(xticksep) "xticksep{" xticksep "} " endif
	if isvar(yticksep) "yticksep{" yticksep "} " endif
	if isvar(xnsubtick) "xnsubtick{" xnsubtick "} " endif
	if isvar(ynsubtick) "ynsubtick{" ynsubtick "} " endif
	if isvar(color) "color{" color "} " endif
	if isvar(xlabel) "belowgap{0.2c} belowcaption{" xlabel "} " endif
	if isvar(ylabel) "leftcaption{ +90d @Rotate{" ylabel "}} " endif
	"objects{"
	return "";
endfunction

function graphcaption(s)
	"} rightcaption {"
	return "";
endfunction

function begingraph()
	return "}{"
endfunction

function endgraph()
	return "}"
endfunction

function graphobj(dir,at,obj)
	"@"|dir " at{" at "} "
	if isvar(margin)
		"margin{" margin "} {" obj "} " 
	else
		"margin{0.5f} {" obj "} "
	endif
	return "";
endfunction

function graph2d(x,y=0)
	if !isvar(color)
		if rows(y)==1 color="red";
		else color=1:rows(y);
		endif
	endif
	if iscomplex(x)
		"@Data "
		if isvar(pairs)
			"pairs{" pairs "} "
		else
			"pairs{solid} "
		endif
		if isvar(linewidth) "linewidth{" linewidth "} " else "linewidth{1.5p} " endif
		if isvar(points)
			"points{" points "} symbollinewidth{1.5p} "
			if isvar(symbolsize) "symbolsize{" symbolsize "} "
			else "symbolsize{0.3f}"
			endif
		endif
		if isstring(color) "color{" color "} "
		else "color{" getcolor(color[i]) "} "
		endif
		"{ "
		re(x)'|im(x)' 
		" }"
		
	else
		for i=1 to rows(y)
			"@Data "
			if isvar(pairs)
				"pairs{" pairs "} "
			else
				"pairs{solid} "
			endif
			if isvar(linewidth) "linewidth{" linewidth "} " else "linewidth{1.5p} " endif
			if isvar(points)
				"points{" points "} symbollinewidth{1.5p} "
				if isvar(symbolsize) "symbolsize{" symbolsize "} "
				else "symbolsize{0.3f}"
				endif
			endif
			if isstring(color) "color{" color "} "
			else "color{" getcolor(color[i]) "} "
			endif
			"{ "
			(x_y[i,:])' 
			" }"
		end
	endif
	return "";
endfunction

function ustep(t)
	return (t>0);
endfunction

function sawwave(t,f)
## returns a vector containing a unit saw wave
##
## t     : time
## f     : frequency
	return t*f-floor(t*f);
endfunction

function sqrwave(t,f,E=1,alpha=0.5)
## returns a vector containing the square signal
##
## t     : time
## f     : frequency
## E     : amplitude (from 0 to max)
## alpha : period ratio
##
	return E*(2*(sawwave(t,f)<alpha)-1);
endfunction;

function triwave(t,f,E=1,alpha=0.5)
## returns a vector containing the triangle signal
##
## t     : time
## f     : frequency
## E     : amplitude (from 0 to max)
## alpha : period ratio
##
	s=(t*f-floor(t*f));
	return E*(2*(s/alpha*(s<=alpha)+(1-s)/(1-alpha)*(s>alpha))-1);
endfunction;

function tf2zp(num,den)
## {z,p,k}=tf2zp(num,den)
## Convert transfer function filter parameter to zero-pole-gain form
## See also : zp2tf, ss2tf, tf2ss, ss2zp, zp2ss.
	return {polysolve(num),polysolve(den),num(length(num))/den(length(den))};
endfunction


function zp2tf(z,p,k)
## {num,den}=zp2tf(z,p,k)
## Convert zero-pole-gain filter parameters to transfer function form.
## See also : tf2zp, ss2tf, tf2ss, ss2zp, zp2ss.
  if length(z)==0 b=1; else b=re(polycons(z)); endif
  return {k*b,re(polycons(p))};
endfunction


function tf2ss(num,den)
## {A,B,C,D} = tf2ss(num,den)
## transforms the transfer function coefficient to state space matrices
## See also : ss2tf, impulse, step.
	k=den[length(den)];
	de=den/k;
	nu=num/k;
	n=length(de)-1;
	if length(nu)>length(de);
		error("the degree of numerator must be less or equal to the denominator's one.");
	endif;
	if length(nu)==length(de);
		{q,r}=polydiv(nu,de);
		D=q;
		nu=r[1:n];
	else
		D=0;
	endif;
	A=zeros([n,n]);
	loop 1 to n-1;
		A[#,#+1]=1;
	end;
	A[n,:]=-de[1:n];
	B=zeros([n,1]);
	B[n,1]=1;
	C=zeros([1,n]);
	C[1:length(nu)]=nu;
	return {A,B,C,D};
endfunction


function impulse()
## {y,x}=impulse(A,B,C,D,t)
## {y,x}=impulse(num,den,t)
## computes the impulse response of the linear system defined in the state
## space (D = 0) or in transfer function form.
## y : response matrix (one row per output)
## x : state variable evolution (one row per state variable)
## For muti-input system, you'll have to consider each input at a time,
## selecting for C one column of the original input matrix each time.
## See also : tf2ss, step.
	if argn()<5; .. assume it is in transfer function form
		{A,B,C,D}=tf2ss(arg1,arg2);t=arg3;
	else;
		A=arg1;B=arg2;C=arg3;D=arg4;t=arg5;
	endif;
	dt=t[2]-t[1];
	s=size(A);n=s[1];s=size(C);m=s[1];
	eAdt=inv(id(n)-A*dt+0.5*(A*dt).(A*dt)); .. Pade approximation
	x=zeros(n,length(t));
	y=zeros(m,length(t));
	x[:,1]=B;
	y[:,1]=C.x[:,1];
	loop 2 to length(t);
		x[:,#]=eAdt.x[:,#-1];
		y[:,#]=C.x[:,#];
	end;
	return {y,x};
endfunction


function step()
## {y,x}=step(A,B,C,D,t,U=1)
## {y,x}=step(num,den,t,U=1)
## computes the unit step response of the linear system defined in the state
## space or in transfer function form.
## U : input step value
## y : response matrix (one row per output)
## x : state variable evolution (one row per state variable)
## For muti-input system, you'll have to consider each input at a time,
## selecting for C one column of the original input matrix each time.
## See also : tf2ss, impulse.
	if argn()<5; .. assume it is in transfer function form
		{A,B,C,D}=tf2ss(arg1,arg2);t=arg3;
		if argn()==3;
			U=1;
		else;
			U=arg4;
		endif;
	else;
		A=arg1;B=arg2;C=arg3;D=arg4;t=arg5;
		if argn()==5;
			U=1;
		else;
			U=arg6;
		endif;
	endif;
	dt=t[2]-t[1];
	s=size(A);n=s[1];s=size(C);m=s[1];
	eAdt=inv(id(n)-A*dt+0.5*(A*dt).(A*dt)); .. Pade approximation
	M=inv(A).(eAdt-id(n)).(B*U);
	x=zeros(n,length(t));
	y=zeros(m,length(t));
	loop 2 to length(t);
		x[:,#]=eAdt.x[:,#-1]+M;
		y[:,#]=C.x[:,#]+D*U;
	end;
	return {y,x};
endfunction


function lsim()
## {y,x}=lsim(A,B,C,D,u,t,[x0])
## {y,x}=lsim(num,den,u,t)
## computes the unit step response of the linear system defined in the state
## space form or in the transfer function form.
## u : input vector
## y : response matrix (one row per output)
## x : state variable evolution (one row per state variable)
## For muti-input system, you'll have to consider each input at a time,
## selecting for C one column of the original input matrix each time.
## See also : tf2ss, impulse, step.
	if argn()<6; .. assume it is in transfer function form
		{A,B,C,D}=tf2ss(arg1,arg2);u=arg3;t=arg4;
		x0=zeros(rows(A),1);
	else;
		A=arg1;B=arg2;C=arg3;D=arg4;u=arg5;t=arg6;
		if argn()==7;
			x0=arg7;
		else;
			x0=zeros(rows(A),1);
		endif;
	endif;
	dt=t[2]-t[1];
	s=size(A);n=s[1];s=size(C);m=s[1];
	eAdt=inv(id(n)-A*dt+0.5*(A*dt).(A*dt)); .. Pade approximation
	M=inv(A).(eAdt-id(n));
	x=zeros(n,length(t));x[:,1]=x0;
	y=zeros(m,length(t));y[:,1]=C.x[:,1]+D*u[:,1];
	loop 2 to length(t);
		x[:,#]=eAdt.x[:,#-1]+M.(B*u[#]);
		y[:,#]=C.x[:,#]+D*u[#];
	end;
	return {y,x};
endfunction

function cheb(x,n)
	signum=-mod(n,2)*2+1;
	w1=matrix([rows(n),0],0);
	w2=matrix([rows(n),0],0);
	w3=matrix([rows(n),0],0);
	if x>1
		i=nonzeros(x>1);
		w1=(x[i]+sqrt(x[i]^2-1))^n;
		w1=(w1+1/w1)/2;
	endif
	if x<-1
		i=nonzeros(x<-1);
		w3=(-x[i]+sqrt(x[i]^2-1))^n;
		w3=signum*(w3+1/w3)/2;
	endif
	if x>=-1 && x<=1
		i=nonzeros(x>=-1 && x<=1);
		w2=cos(n*acos(x[i]));
	endif
	return w3|w2|w1;
endfunction

function alp1(f,fo)
## Passe bas asymptotique 1er ordre (asymptotical low pass)
## fo : fréquence de cassure
## f : vecteur fréquence
	g=ones(size(f));
	.. Calcul des fréquences après la fréquence de coupure
	i=nonzeros(f>fo);
	g[i]=-1i*fo/f[i];
	return g;
endfunction

function ahp1(f,fo)
## Passe haut asymptotique 1er ordre (asymptotical high pass)
## fo : fréquence de cassure
## f  : vecteur fréquence
	g=ones(size(f));
	.. Calcul des fréquences après la fréquence de coupure
	i=nonzeros(f<fo);
	g[i]=1i*f[i]/fo;
	return g;
endfunction

function alp2(f,fo)
## Passe bas asymptotique 2e ordre (asymptotical low pass)
## fo : fréquence de cassure
## f  : vecteur fréquence
	g=ones(size(f));
	.. Calcul des valeurs après la fréquence de coupure
	i=nonzeros(f>=fo);
	g[i]=-(fo/f[i])^2;
	return g;
endfunction

function ahp2(f,fo)
## Passe haut asymptotique 2e ordre (asymptotical high pass)
## fo : fréquence de cassure
## f  : vecteur fréquence
	g=ones(size(f));
	.. Calcul des valeurs avant la fréquence de coupure
	i=nonzeros(f<fo);
	g[i]=-(f[i]/fo)^2;
	return g;
endfunction


function butterord()

endfunction

function butter(n,wc)
## butter(n,wc)
##  returns the Butterworth transfert function in zpk form
  return {[],wc*exp(-1i*pi*(2*(0:n-1)+n+1)/(2*n)),wc^n};
endfunction


function cheb1(n,wc,rp)
  eps=sqrt(10^(rp/10)-1);
  g=((1+sqrt(1+eps^2))/eps)^(1/n);
  i=1:n;
  p=wc*((1/g-g)/2*sin((2*i-1)*pi/(2*n))+1i*(1/g+g)/2*cos((2*i-1)*pi/(2*n)));
  k=re(prod(-p));
  if mod(n,2)==0 k=k*10^(-rp/20); endif
  return {[],p,k};
endfunction


function ellipord(wp,ws,Ap,As)
  k=wp/ws;  .. Selectivity factor k
  u=(1-(1-k^2)^0.25)/(2*(1+(1-k^2)^0.25));
  q=u+2*u^5+15*u^9+150*u^13;       .. modular constant q
  D=(10^(As/10)-1)/(10^(Ap/10)-1); .. discrimimnation factor D
  n=ceil(log(16*D)/log(1/q));      .. minimum filter order
  Ts=10*log10(1+(10^(Ap/10)-1)/(16*q^n)); .. effective minimum stopband loss
  return {n,Ts}
endfunction

function ellip(n,wp,ws,Ap)
  k=wp/ws;alpha=sqrt(wp*ws);
  u=(1-(1-k^2)^0.25)/(2*(1+(1-k^2)^0.25));
  q=u+2*u^5+15*u^9+150*u^13;       .. modular constant q
  N=4;
  V=1/(2*n)*log((10^(Ap/20)+1)/(10^(Ap/20)-1));
  m=0:N;l=1:N;
  p0=abs(q^0.25*sum((-1)^m*q^(m*(m+1))*sinh((2*m+1)*V))/(0.5+sum((-1)^l*q^(l^2)*cosh(2*l*V))));
  W=sqrt((1+p0^2/k)*(1+k*p0^2));
  r=floor(n/2);
  t=(2*r~=n);
  mu=((1:r)*(1-t)+((1:r)-0.5)*t)';
  m=0:N;l=1:N;
  X=2*q^0.25*sum((-1)^m*q^(m*(m+1))*sin((2*m+1)*mu*pi/n))/(1+2*sum((-1)^l*q^(l^2)*cos(2*l*mu*pi/n)));
  Y=sqrt((1-X^2/k)*(1-k*X^2));
  a=1/X^2;
  b=2*p0*Y/(1+(p0*X)^2);
  c=((p0*Y)^2+(X*W)^2)/(1+(p0*X)^2)^2;
  H0=prod((c')/(a'))*(10^(-Ap/20)*t+p0*(1-t));
  den=c|b|ones(size(b));
  p=[];z=[];
  if t~=0 p=p|(-p0); endif
  for i=1 to r
    p=p|polysolve(den[i,:]);
    z=z|[-1i*sqrt(a[i]), 1i*sqrt(a[i])];
  end
  z=z*alpha;p=p*alpha;H0=H0*(alpha*(length(p)-length(z))*(1-t)+t);
  return {z,p,H0};
endfunction

function be(n)
  if n==1 return [1,1]; endif
  qn0=1;qn1=[1,1,0];
  loop 2 to n;
    qn=(2*#-1)*qn1+([0,0]|qn0);
    qn0=qn1[1:#];
    qn1=qn|0;
  end;
  return qn;
endfunction

function renorm(b)
  a=[1,0.72675,0.57145,0.46946,0.41322,0.37038,0.33898,0.31546];
  n=length(b)-1;
  if n<=8;
    loop 1 to n;
      b[#]=a[n]^(n-#+1)*b[#];
    end;
  endif;
  return b;  
endfunction

function bessel(n,wc,norm="delay")
  normfactor=[1,0.72675,0.57145,0.46946,0.41322,0.37038,0.33898,0.31546];
  b=be(n);k=0;
  p=polysolve(b);p=p*wc;
  k=b[1]*wc^n;
  if norm=="mag"
    p=p*normfactor[n];
    k=k*(normfactor[n])^n;
  endif
  return {[],p,k};
endfunction

function reverse(a)
## reverse(p)
##   reverse a vector
  v=a;
  loop 1 to length(a);
    v[#]=a[length(a)-#+1];
  end;
  return v;
endfunction

function freqs()
## computes the frequency response (w=pulsation)
## H=freqs(b,a,w) 
##   transfer function in the form b/a
## H=freqs(z,p,k,w)
##   transfer function in the form zeros/poles/gain
  if argn==3
    return polyval(arg1,1i*arg3)/polyval(arg2,1i*arg3);
  elseif argn()==4
  	if length(arg1)==0
  	  num=1;
  	else
      num=re(polycons(arg1));
    endif
    den=re(polycons(arg2));
    return arg3*polyval(num,1i*arg4)/polyval(den,1i*arg4);
  endif
  error("bad parameters");
  return 0;
endfunction

function mag()
## computes the frequency response magnitude (w=pulsation)
## H=freqs(b,a,w) 
##   transfer function in the form b/a
## H=freqs(z,p,k,w)
##   transfer function in the form zeros/poles/gain
  if argn==3
    return abs(polyval(arg1,1i*arg3)/polyval(arg2,1i*arg3));
  elseif argn()==4
  	if length(arg1)==0
  	  num=1;
  	else
      num=re(polycons(arg1));
    endif
    den=re(polycons(arg2));
    return abs(arg3*polyval(num,1i*arg4)/polyval(den,1i*arg4));
  endif
  error("bad parameters");
  return 0;
endfunction

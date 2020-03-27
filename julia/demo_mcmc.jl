#simple MCMC to sample a Gaussian in Julia
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using UseGA_AEM, Random

#forward

#system geometry constants
const ztxLM = 35.0;
const ztxHM = 35.0;
const rrx = -17.0;
const dzrxLM = 2.0;
const dzrxHM = 0.2;

#cpp struct geometry
g = UseGA_AEM.Geometry(ztxLM = ztxLM,
						   ztxHM = ztxHM,
						   rrx = rrx,
						   dzrxLM = dzrxLM,
						   dzrxHM = dzrxHM);

#initialise AEM system (other params are hard-wired in the
#C++ or the julia module for now)
em = UseGA_AEM.init_GA_AEM();

#forward function
f = UseGA_AEM.EMoperator(em,g);
function forward(x)
	f(ztxLM,ztxHM,x[1:length(x)÷2+1],x[length(x)÷2+1:end])
	[f.em.SZLM;f.em.SZHM]
end

#make synth data
function noisy_forward(x; nmag=3.0)
	signal = forward(x)
	noise = 1.0 .+ nmag/100 .* randn(size(signal))
	signal .* noise
end

#functions
function MCMCstep(x,misfit,params)
	priorViolate = false;
	accept = false;
	xNew = copy(x);
	k = rand(1:length(x));
	xNew[k] += params["rSD"][k]*randn();
	if xNew[k] < params["rmin"][k] || xNew[k] > params["rmax"][k]
		priorViolate = true;
	end

	if !priorViolate
		newmisfit = get_misfit(xNew);
		logalpha = -1/params["T"] * (newmisfit - misfit);
		# println(logalpha)
		if log(rand()) < logalpha
			x = xNew;
			accept = true;
			misfit = newmisfit;
		end
	end

	x, misfit, accept
end

Random.seed!(2);
true_c = [0.01,1.0,0.01];
true_t = [50,50];
noisy_data = noisy_forward([true_c;true_t]);
σ2 = (0.03*abs.(noisy_data)).^2;

function get_misfit(logc_x,nmag=σ2)
	x = copy(logc_x);
	x[1:length(x)÷2+1] = 10 .^ logc_x[1:length(logc_x)÷2+1];
	response = forward(x);
	res = response - noisy_data;
	sum((res.^2)./(2*σ2))
end

function runMCMC(nsamples)

	nlayers = 3;

	#parameters for the proposal and prior
	params = Dict()

	#width of proposal distribution for conductivity and thickness
	params["rSD"] = [0.1*ones(nlayers);10*ones(nlayers-1)];

	#prior is a uniform dist over some bounds. log for conductivity
	params["rmin"] = [-3*ones(nlayers);zeros(nlayers-1)];
	params["rmax"] = [1.5*ones(nlayers);100*ones(nlayers-1)];


	Random.seed!(2);

	#initialise the chain state by sampling from the prior
	x = params["rmin"] + (params["rmax"]-params["rmin"]).*rand(length(params["rSD"]),1);
	#can change the temperature to implement PT if desired
	params["T"] = 1.0;
	#track acceptance ratio and aim for 20-ish percent
	accepted = 0;
	misfit = get_misfit(x);

	#store samples and misfits
	max_depth = 300;
	ncells = 61;
	X = zeros(nsamples,ncells);
	Xct = zeros(nsamples,length(x))
	chi2by2 = zeros(nsamples,1);

	#run the loop. This obviously must be serial because next sample
	#depends on current state.
	for i = 1:nsamples
		x, misfit, accept = MCMCstep(x,misfit,params);
		#store the model as conductivity-depth for easier analysis
		X[i,:] = ct2cd(x,max_depth,ncells);
		Xct[i,:] = x;
		chi2by2[i] = misfit;
		if accept
			accepted += 1;
		end

	end

	println(accepted/nsamples)

	X,Xct,chi2by2
end

function ct2cd(x,max_depth,ncells)
	c = x[1:length(x)÷2 + 1];
	#interface depths
	st = cumsum(x[length(x)÷2 + 2:end]);

	#fill in each layer
	cdarr = c[end]*ones(ncells);
	d = range(0,max_depth,length=ncells);
	for i = length(st):-1:1
		cdarr[d .< st[i]] .= c[i]
	end

	cdarr
end
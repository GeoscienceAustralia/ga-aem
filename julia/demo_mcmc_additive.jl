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
const max_depth = 150.0;

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
	f(ztxLM,ztxHM,x[1:length(x)÷2+1],x[length(x)÷2+2:end]);
	[f.em.SZLM;f.em.SZHM]
end

#make synth data with additive noise
function noisy_forward(x,anmag; nmag=0.03)
	signal = forward(x);
	noise = 1.0 .+ nmag .* randn(size(signal));

	#don't want it being negative
	abs.((signal .* noise) .+ anmag .* randn(size(signal)));
end

Random.seed!(2);
true_c = [0.01,1.0,0.01];
true_t = [50,50];

true_response = forward([true_c;true_t]);
ndata = length(true_response);

#additive noise magnitudes
amag = 1e-13 * ones(ndata);

noisy_data = noisy_forward([true_c;true_t],amag);

rms_err = sqrt(sum(((noisy_data-true_response)./noisy_data).^2)/ndata);
println(rms_err);

#put the additive in same units as multiplicative for convenience.
avar = (amag ./ noisy_data).^2;

#this is the maximally-costly misfit calculator.
#it does a full computation of the forward model to compute residuals.
#use the noisemove_misfit_ratio function for noise magnitude changes 
#(or noisemove_with_additive if you want additive noise as well)
#noise variance is treated as a vector here so this works for additive noise as well
function get_misfit(logc_x,var,d2vec=noisy_data.^2)
	x = copy(logc_x);
	x[1:length(x)÷2+1] = 10 .^ logc_x[1:length(logc_x)÷2+1];
	response = forward(x);
	res = response - noisy_data;
	sum((res.^2)./(2*var .* d2vec)), (res.^2)./d2vec
end

#functions
function MCMCstep(x,nmag,var,misfit,res,params;max_depth=max_depth)
	priorViolate = false;
	accept = false;
	xNew = copy(x);
	k = rand(1:length(x)+1);
	#'normal' moves
	if k <= length(x)
		xNew[k] += params["rSD"][k]*randn();
		if xNew[k] < params["rmin"][k] || xNew[k] > params["rmax"][k]
			priorViolate = true;
		end

		if sum(x[length(x)÷2+2:end])>max_depth
			priorViolate = true;
		end

		if !priorViolate
			newmisfit, newres = get_misfit(xNew,var);
			logalpha = -1/params["T"] * (newmisfit - misfit);
			# println(logalpha)
			if log(rand()) < logalpha
				x = xNew;
				accept = true;
				misfit, res = newmisfit, newres;
			end
		end
	else
		nmag_new = nmag + params["nSD"]*randn();
		if nmag_new < params["nmin"] || nmag_new > params["nmax"]
			priorViolate = true;
		end

		if !priorViolate
			var_new = var .- nmag^2 .+ nmag_new^2;
			logalpha = -1/params["T"] * noisemove_with_additive(res,var,var_new);
			if log(rand()) < logalpha
				#chi^2 scales quite simply - don't need the normalising term for this
				misfit = sum(res./(2*var_new));
				nmag = nmag_new;
				var = var_new;
				accept = true;
			end
		end

	end

	x, nmag, var, misfit, res, accept
end

function runMCMC(nsamples;max_depth=max_depth)

	nlayers = 3;

	#parameters for the proposal and prior
	params = Dict()

	#width of proposal distribution for conductivity and thickness
	params["rSD"] = [0.5*ones(nlayers);25*ones(nlayers-1)];

	#prior is a uniform dist over some bounds. log for conductivity
	params["rmin"] = [-3*ones(nlayers);25*ones(nlayers-1)];
	params["rmax"] = [1*ones(nlayers);100*ones(nlayers-1)];

	#initialise min and max for noise variance proposal
	params["nmin"] = 0.01;
	params["nmax"] = 0.1;
	params["nSD"] = 0.02;

	Random.seed!();

	#initialise the chain state by sampling from the prior
	tooDeep = true;
	x = zeros(2*nlayers-1);
	while tooDeep
		x = params["rmin"] + (params["rmax"]-params["rmin"]).*rand(length(params["rSD"]),1);
		if sum(x[nlayers+1:end]) <= max_depth
			tooDeep = false;
			println(x)
		end
	end

	#initialise noise magnitude estimate
	nmag = params["nmin"] + (params["nmax"]-params["nmin"])*rand();
	var = avar .+ nmag^2;

	#can change the temperature to implement PT if desired
	params["T"] = 1.0;
	#track acceptance ratio and aim for 20-ish percent
	accepted = 0;
	misfit, res = get_misfit(x,var);

	#store samples and misfits
	max_depth = 300;
	ncells = 61;
	X = zeros(nsamples,ncells);
	Xct = zeros(nsamples,length(x))
	chi2by2 = zeros(nsamples,1);
	nmag_vec = zeros(nsamples,1);

	#run the loop. This obviously must be serial because next sample
	#depends on current state.
	for i = 1:nsamples
		x, nmag, var, misfit, res, accept = MCMCstep(x,nmag,var,misfit,res,params,max_depth=max_depth);
		#store the model as conductivity-depth for easier analysis
		X[i,:] = ct2cd(x,max_depth,ncells);
		Xct[i,:] = x;
		chi2by2[i] = misfit;
		nmag_vec[i] = nmag;
		if accept
			accepted += 1;
		end

	end

	println(accepted/nsamples)

	X,Xct,chi2by2, nmag_vec
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

#this needs to be a separate function because changing noise amplitude
#adds a normalising term to the ratio of likelihoods, and also doesn't
#require re-computing the forward
function noisemove_misfit_ratio(misfit,nmag,nmag_new,ndata)
	ndata*(log(nmag_new)-log(nmag)) + misfit*((nmag/nmag_new)^2 - 1)
end

function noisemove_with_additive(residuals,var,var_new)
	1/2 * sum(log.(var_new) - log.(var) +  (1 ./ var_new - 1 ./ var).*residuals)
end
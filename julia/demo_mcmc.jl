#simple MCMC to sample a Gaussian in Julia
using Random


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
		newmisfit = get_misfit(xNew,params);
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

function get_misfit(x,params)
	res = x - params["mu"];
	m = (res' * params["Ci"] * res)/2;
	m[1]
end

function runMCMC()
	nsamples = 100000;

	#parameters of the distribution to be sampled
	params = Dict()
	params["C"] = [1 0.6;0.6 0.5];
	params["Ci"] = inv(params["C"]);
	params["mu"] = [1;1];

	#width of proposal distribution

	params["rSD"] = [1;1];

	#prior is a uniform dist over some bounds
	params["rmin"] = [-2;-2];
	params["rmax"] = [4;4];

	Random.seed!(2);

	#initialise the chain state by sampling from the prior
	x = params["rmin"] + (params["rmax"]-params["rmin"]).*rand(length(params["rSD"]),1);
	#can change the temperature to implement PT if desired
	params["T"] = 1.0;
	#track acceptance ratio and aim for 20-ish percent
	accepted = 0;
	misfit = get_misfit(x,params);

	#store samples and misfits
	X = zeros(nsamples,length(x));
	chi2by2 = zeros(nsamples,1);

	#run the loop. This obviously must be serial because next sample
	#depends on current state.

	for i = 1:nsamples
		x, misfit, accept = MCMCstep(x,misfit,params);
		X[i,:] = x;
		chi2by2[i] = misfit;
		if accept
			accepted += 1;
		end

	end

	println(accepted)

	X,chi2by2
end



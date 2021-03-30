with(RegularChains):
with(SemiAlgebraicSetTools):
with(ListTools):

EvalCertification := module()
option package;

export HolderInformationForPolynomial, 
	   HolderInformationForExponential, 
	   HolderInformationForRationalPolynomial, 
	   EstimateRootsAndCertifyEvaluations,
	   ExpandToImaginaryAndReal;

#ExpandToImaginaryAndReal,
local SqrFreeFactorization,
	  SingleEqualityToEstimate,
	  BoxEqualitiesToEstimate,
	  EstimateAndClassifySolutions,
	  HolderDataToCertificationTolerance,
	  EvaluateFunctionAtRoot;

# Compute a square free factorization of a polynomial 
# along with multiplicities
SqrFreeFactorization := proc(F) 
	local FsqrFree, fa, CurrOutput, j, k; FsqrFree := sqrfree(F); 
	FsqrFree := FsqrFree[2]; CurrOutput := []; 
	for j to nops(FsqrFree) do 
		fa := ifelse(type(factor(FsqrFree[j][1]), `+`), [FsqrFree[j][1]], convert(factor(FsqrFree[j][1]), list)); 
		for k to nops(fa) do 
			CurrOutput := [op(CurrOutput), [fa[k], FsqrFree[j][2]]]; 
		end do; 
	end do; 
	return CurrOutput; 
end proc;

# Expand a list of polynomials to list with real and 
# imaginary parts separated
ExpandToImaginaryAndReal := proc(F) 
	local vars, imaginary_variables, expanded_equations, i; 
	vars := indets(F); 
	imaginary_variables := [seq(b[i], i = 1 .. nops(vars))]; 
	expanded_equations := expand(subs([seq(vars[i] = Complex(vars[i], imaginary_variables[i]), i = 1 .. nops(vars))], F)); 
	return [[op(subs(I = 0, expanded_equations)), op([seq(coeff(expanded_equations[i], I), i = 1 .. nops(expanded_equations))])], [seq(vars[i]^2 + imaginary_variables[i]^2 - 1, i = 1 .. nops(vars))]]; 
end proc;

SingleEqualityToEstimate := proc(equality) 
	return ifelse(type(rhs(equality), 'list'), 1/2*rhs(equality)[2] + 1/2*rhs(equality)[1], rhs(equality)); 
end proc;

BoxEqualitiesToEstimate := proc(equalities) 
	local reduced_equalities; 
	reduced_equalities := map(SingleEqualityToEstimate, equalities); 
	return Complex(reduced_equalities[1], reduced_equalities[2]); 
end proc;


# Estimate the solution of univarite polynomial F and classify 
# by whether the imaginary norm is > 1, < 1, or = 1
EstimateAndClassifySolutions := proc(F, epsilon) 
	# INPUT: Squarefree polynomial F in format [polynomial,power], 
	# estimation error epsilon

	# OUTPUT: [list_of_estimates_for_norm_1_solutions,list_of_estimates_for_norm_<_1_solutions,
	# 		   list_of_estimates_for_norm_>_1_solutions]
	local actual_polynomial, system_with_imaginary_parts, equations, 
		  norms, R, norm_1_boxes, norm_1_values, norm_greater_boxes, 
		  norm_greater_values, norm_less_boxes, norm_less_values, i, 
		  norm_1_values_refined, norm_greater_values_refined, norm_less_values_refined; 
	actual_polynomial := [F[1]]; 
	system_with_imaginary_parts := ExpandToImaginaryAndReal(actual_polynomial); 
	equations := system_with_imaginary_parts[1]; 
	norms := system_with_imaginary_parts[2]; 
	R := PolynomialRing([seq(indets(norms)[i], i = 1 .. nops(indets(norms)))]); 
	norm_1_boxes := RealRootIsolate([op(equations)], [], [], [], R, 'abserr' = epsilon); 
	norm_1_values := map(BoxValues, norm_1_boxes, R); 
	norm_1_values_refined := map(BoxEqualitiesToEstimate, norm_1_values); 
	norm_1_values_refined := [seq(Record("value"=norm_1_values_refined[i]), i = 1 .. nops(norm_1_values_refined))]; 
	return [op(norm_1_values_refined)]; 
end proc;

# Compute local Lipschitz constant for univarite polynomial F at a point
# given an estimate for that point and a bound on that estimate's accuracy
HolderInformationForPolynomial := proc(F, PointEstimate, EstimateAccuracy)
	# INPUT: Univariate polynomial F, 
	# an estimate PointEstimate of a point in F's domain (could be real or complex),
	# an upper bound EstimateAccuracy on the distance from PointEstimate to the
	# actual root

	# OUTPUT: A Record of the form [exponent,constant] such that
	# F is certifiably Holder with those values
	local deg, Derivatives, HolderConstant, HolderExponent, i; 
	deg := degree(F); 
	if deg > 0 then		
		Derivatives := Vector();
		Derivatives(1) := diff(F, indets(F)[1] $ 1);
		for i from 2 to deg do
	    	Derivatives(i) := diff(Derivatives(i - 1), indets(F)[1] $ 1);
		end do;
		Derivatives := convert(Derivatives, list); 
		HolderConstant := add((EstimateAccuracy)^(i - 1)*abs(subs(indets(F)[1] = PointEstimate, Derivatives[i]))/(i - 1)!, i = 1 .. deg);
	else HolderConstant := 0;
	end if;
	HolderExponent := 1; 
	return Record('exponent'=HolderExponent,'constant'=HolderConstant,'avoid_roots'=false,'max_degree'=deg-1);
end proc;

# Holder information for exponentials of the form x^alpha
HolderInformationForExponential := proc(function_exponent) 
	local output_function; 
	output_function := proc(F, PointEstimate, EstimateAccuracy)
		return Record('exponent'=function_exponent,'constant'=1,'avoid_roots'=false,'max_degree'=1);
	end proc; 
	return output_function; 
end proc;

HolderInformationForRationalPolynomial := proc(P, PointEstimate, EstimateAccuracy,domain_estimate:=false) 
    # INPUT: Univariate rational polynomial  of the form F/G

    # OUTPUT: A record of the form [exponent, constant] such that 
    # F/G is certifiably Holder with those values
	local F, G, differential_numerator, numerator_constant, max_value, 
		  denominator_constant, min_value,FactorsAndMultiplicities,ClassifiedRoots,MINIMUM_NONZERO; 
	if domain_estimate then		
		FactorsAndMultiplicities := SqrFreeFactorization(denom(simplify(P)));
		ClassifiedRoots := map(EstimateAndClassifySolutions, FactorsAndMultiplicities, EstimateAccuracy);
		ClassifiedRoots := FlattenOnce(ClassifiedRoots);
		return Record('avoid_roots'=[seq(ClassifiedRoots[i]:-value,i=1..nops(ClassifiedRoots))]);
	end if;

	MINIMUM_NONZERO := 3*10^(-10);
	F := numer(simplify(P)); 
	G := denom(simplify(P)); 
	
	if type(F, integer) then 
		differential_numerator := -diff(G, indets(G)[1] $ 1)*F; 
	elif type(G, integer) then 
		differential_numerator := G*diff(F, indets(F)[1] $ 1); 
	else differential_numerator := -diff(G, indets(G)[1] $ 1)*F + G*diff(F, indets(F)[1] $ 1); 
	end if; 	
	numerator_constant := HolderInformationForPolynomial(differential_numerator, PointEstimate, EstimateAccuracy):-constant; 

	if type(differential_numerator, integer) then 
		max_value := abs(differential_numerator) + numerator_constant*EstimateAccuracy; 
	else max_value := abs(subs(indets(differential_numerator)[1] = PointEstimate, differential_numerator)) + numerator_constant*EstimateAccuracy; 
	end if; 
	if type(G, integer) then 
		denominator_constant := 1; 
	else denominator_constant := HolderInformationForPolynomial(G^2, PointEstimate, EstimateAccuracy):-constant; 
	end if; 
	if evalb(numerator_constant = denominator_constant) then
		return MINIMUM_NONZERO;
	end if;
	if type(G, integer) then
		min_value := G^2 - denominator_constant*EstimateAccuracy; 
	else min_value := abs(subs(indets(G)[1] = PointEstimate, G))^2 - denominator_constant*EstimateAccuracy; 
	end if; 

	return Record('exponent'=1,'constant'=max_value/min_value,'max_degree'=1);
end proc;

HolderDataToCertificationTolerance := proc(HolderData, ErrorTolerance) 
	local HolderConstant, HolderExponent,RoundedDigitsTolerance,RoundedCertNegativeExponent;
	HolderExponent := HolderData:-exponent;
	HolderConstant := HolderData:-constant;
	RoundedDigitsTolerance := max(-floor(log[10](ErrorTolerance^HolderData:-max_degree)),1) + 1;
	if evalf[RoundedDigitsTolerance](HolderConstant) <= 2*10^(-RoundedDigitsTolerance) then
		return 0;
	end if;
	RoundedCertNegativeExponent := floor(-log[10](ErrorTolerance/HolderConstant)/HolderExponent);
	return min(1/10^RoundedCertNegativeExponent, ErrorTolerance, 1);
end proc;

EvaluateFunctionAtRoot := proc(Function, Root, Precision)
	local SymbolicEval,DigitsBeforeDecimalPoint,DigitsAfterDecimalPoint;
	if type(Function,procedure) then
		SymbolicEval := Function(Root);
	else
		SymbolicEval := subs(indets(Function)[1]=Root,Function);
	end if;
	if evalb(SymbolicEval=0) then
		return 0;
	end if;
	# Precision calculation: Figure out necessary number of digits
	# in front of the decimal point, figure out the same for after
	# the decimal point, then add them together.
	DigitsBeforeDecimalPoint := max(ceil(evalf[1](log[10](abs(SymbolicEval)))),1); 
	DigitsAfterDecimalPoint := max(-floor(log[10](Precision)),0);
	return evalf[max(1,DigitsBeforeDecimalPoint+DigitsAfterDecimalPoint+1)](SymbolicEval);
end proc;

EstimateRootsAndCertifyEvaluations := proc(PolynomialToSolve, FunctionsToEvaluate, HolderEstimationProcedure, ErrorTolerance, AbandonThreshold := 10^(-100)) 
	local FactorsAndMultiplicities, Multiplicities, ClassifiedRoots, 
	      RootsWithoutAdditionalInformation, HolderBounds, HolderCertificationBounds, 
	      SharpenedRootTolerance, ClassifiedRootsWithMultiplicities, WorkingTolerance, ActualTolerance, 
	      i, j,arguments_for_output,root_information_for_output,evaluation_information_for_output,current_evals,
	      minimum_gap,current_roots,holder_function,function_roots,estimated_gap_values;
	FactorsAndMultiplicities := SqrFreeFactorization(expand(PolynomialToSolve));
	Multiplicities := [seq(FactorsAndMultiplicities[i][2], i = 1 .. nops(FactorsAndMultiplicities))];
	SharpenedRootTolerance := 0;
	# We need to make sure all estimates are at least 1/10th sharper 
	# than the desired number because rounding error in the last reported digit
	# could result in that amount of error
	ActualTolerance := ErrorTolerance - (1/10)*ErrorTolerance;
	WorkingTolerance := max(1/10,ActualTolerance);

	# First Holder information step: Tightening the precision
	# to avoid places where the evaluation functions aren't 
	# defined. We need that there are no undefined points
	# with 2*WorkingTolerance of the function to solve's
	# estimated roots 
	for i to nops(FunctionsToEvaluate) do
		minimum_gap := 0;
		if type(HolderEstimationProcedure,list) then 
			holder_function := HolderEstimationProcedure[i];
		else
			holder_function := HolderEstimationProcedure;
		end if;
		# Note: Minimum gap often contains square roots and Maple doesn't 
		# automatically compare expressions with abs in them symbolically.
		# For minimum_gap we can just square, for others we need to evaluate 
		# and be careful with precision
		while minimum_gap^2 <= (6*WorkingTolerance)^2 and WorkingTolerance > AbandonThreshold do
			if minimum_gap^2 > 0 then
				WorkingTolerance := min(WorkingTolerance/2,minimum_gap/(6*WorkingTolerance));
			end if;
			function_roots := holder_function(FunctionsToEvaluate[i],1.1,WorkingTolerance,true):-avoid_roots;	
			if evalb(function_roots=false) then 
				break;
			end if;	
			ClassifiedRoots := map(EstimateAndClassifySolutions, FactorsAndMultiplicities, WorkingTolerance);
			ClassifiedRoots := FlattenOnce(ClassifiedRoots);
			current_roots := [seq(ClassifiedRoots[j]:-value, j = 1 .. nops(ClassifiedRoots))];
			estimated_gap_values := FlattenOnce([ seq([seq(abs(function_roots[j] - current_roots[k]),k=1..nops(current_roots))],j=1..nops(function_roots)) ]);
			minimum_gap := min(estimated_gap_values);
			# Need to round down to some rational representation 
			# to keep things rational
			minimum_gap := convert(evalf[10](minimum_gap)-2*10^(-10),rational);			
		end do 
	end do;

	# At this point we've reduced to a good enough tolerance to 
	# find local Holder information while dodging points not in the domains 
	# of the functions to evaluate. Proceed to find that information, sharpen tolerance, 
	# and give final results.
	while SharpenedRootTolerance = 0 and AbandonThreshold < WorkingTolerance do 
		ClassifiedRoots := map(EstimateAndClassifySolutions, FactorsAndMultiplicities, WorkingTolerance);
		ClassifiedRoots := FlattenOnce(ClassifiedRoots);
		RootsWithoutAdditionalInformation := [seq(ClassifiedRoots[i]:-value, i = 1 .. nops(ClassifiedRoots))];
		if type(HolderEstimationProcedure, list) then
			if nops(HolderEstimationProcedure) <> nops(FunctionsToEvaluate) then 
				error "Number of Holder estimation procedures is different than the number of functions to evaluate.";
			end if;
			HolderBounds := [seq(map[2](HolderEstimationProcedure[i], FunctionsToEvaluate[i], RootsWithoutAdditionalInformation, 2*WorkingTolerance), i = 1 .. nops(FunctionsToEvaluate))];
		else
			HolderBounds := [seq(map[2](HolderEstimationProcedure, FunctionsToEvaluate[i], RootsWithoutAdditionalInformation, 2*WorkingTolerance), i = 1 .. nops(FunctionsToEvaluate))];
		end if;
		HolderBounds := FlattenOnce(HolderBounds);
		HolderCertificationBounds := map(HolderDataToCertificationTolerance, HolderBounds, min(ActualTolerance,WorkingTolerance/2));
		SharpenedRootTolerance := min(HolderCertificationBounds);
		WorkingTolerance := 1/10*WorkingTolerance;
	end do;
	if SharpenedRootTolerance = 0 then 
		error "Insufficient precision available to certify evaluations.";
	end if;
	ClassifiedRoots := map(EstimateAndClassifySolutions, FactorsAndMultiplicities, SharpenedRootTolerance);
	ClassifiedRoots := FlattenOnce(ClassifiedRoots);
	arguments_for_output := Record("polynomial"=PolynomialToSolve,"evaluated_functions"=FunctionsToEvaluate,"error_tolerance"=ActualTolerance);
	root_information_for_output := Record("root_values"=[seq(ClassifiedRoots[i]:-value, i=1..nops(ClassifiedRoots))], "root_multiplicities"=Multiplicities);
	evaluation_information_for_output := Record();
	for i to nops(FunctionsToEvaluate) do
		current_evals := Record("evaluations_functions_"||i=map[2](EvaluateFunctionAtRoot,FunctionsToEvaluate[i],root_information_for_output:-root_values,ActualTolerance));
		evaluation_information_for_output := Record[evaluation_information_for_output,current_evals]();
	end do;
	return Record[arguments_for_output,root_information_for_output,evaluation_information_for_output]();
end proc;

end module:
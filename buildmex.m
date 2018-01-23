N = 32768;
Nrecur = 32;

cd utils

if ~verLessThan('matlab', '9.0')
	% Uses recursion, so only works in never versions of matlab
	codegen -o genschur genschur.m -args {coder.typeof(complex(1,1), [1 N], [false true]), coder.typeof(complex(1,1), [1 N], [false true])}
else
	% As a fallback we can use mex only for the innermost iteration of the
	% generalized Schur algorithm. This will be slower than the above.
	codegen -o genschur_recurrence genschur_recurrence.m -args {coder.typeof(complex(1,1), [1 Nrecur], [false true]), coder.typeof(complex(1,1), [1 Nrecur], [false true])}
end

cd ..


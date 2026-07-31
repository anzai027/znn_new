function problem = example1_problem()

I2 = eye(2);
Q0 = [I2;-I2];

problem.G  = @(t) [0.25*sin(t) + 1, 0.5*cos(t); ...
                   0.5*cos(t),       0.25*sin(t) + 1];
problem.dG = @(t) [0.25*cos(t), -0.5*sin(t); ...
                  -0.5*sin(t),  0.25*cos(t)];

problem.h  = @(t) [sin(3*t); cos(3*t)];
problem.dh = @(t) [3*cos(3*t); -3*sin(3*t)];

% sin(4t)x1 - cos(4t)x2 = cos(2t)
problem.P  = @(t) [sin(4*t), -cos(4*t)];
problem.dP = @(t) [4*cos(4*t), 4*sin(4*t)];
problem.u  = @(t) cos(2*t);
problem.du = @(t) -2*sin(2*t);

problem.Q  = @(t) Q0; 
problem.dQ = @(t) zeros(4,2); 
problem.v  = @(t) 1e8*ones(4,1); 
problem.dv = @(t) zeros(4,1); 

end

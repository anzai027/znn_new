function problem = problemof(job)
% 按题目设置保留或删除不存在的约束。
if job.example == 2
    problem = example2_problem();
else
    problem = example1_problem();
end
if job.variant == "reduced"
    problem.Q = @(t) zeros(0,2);
    problem.dQ = @(t) zeros(0,2);
    problem.v = @(t) zeros(0,1);
    problem.dv = @(t) zeros(0,1);
end
end

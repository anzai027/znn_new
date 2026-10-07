function list = series(base)
% 固定主尺度并排列补充试验。
list = repmat(base(10),0,1);
for eps = [1e-4,1e-3,1e-2,1e-1]
    for kind = ["zero","cosine"]
        for model = ["VFCR-ZNN","GNN"]
            j = base(find([base.group] == "main" & [base.noise] == kind ...
                & [base.model] == model & [base.seed] == 1,1));
            push(j,eps,"scan");
        end
    end
end
for model = ["VFCR-ZNN","GNN"]
    j = base(find([base.group] == "initial" & [base.scale] == 1 ...
        & [base.model] == model,1));
    push(j,1e-3,"initial");
    j = base(find([base.group] == "precision" & [base.model] == model,1));
    push(j,1e-3,"precision");
end
for seed = 1:5
    for model = ["VFCR-ZNN","GNN"]
        j = base(find([base.group] == "main" & [base.noise] == "gaussian" ...
            & [base.model] == model & [base.seed] == seed,1));
        push(j,1e-3,"random");
    end
end
for rep = 2:3
    for kind = ["zero","cosine"]
        for model = ["VFCR-ZNN","GNN"]
            j = base(find([base.group] == "timing" & [base.noise] == kind ...
                & [base.model] == model & [base.repeat] == rep,1));
            push(j,1e-3,"timing");
        end
    end
end
ids = find([list.eps] == 1e-3 & [list.group] ~= "scan");
for id = ids
    j = list(id);
    j.id = j.base;
    push(j,1e-2,j.group);
end

    function push(j,eps,group)
        j.base = j.id;
        j.id = numel(list)+1;
        j.eps = eps;
        j.group = group;
        if isempty(list)
            list = j;
        else
            list(end+1,1) = j;
        end
    end
end

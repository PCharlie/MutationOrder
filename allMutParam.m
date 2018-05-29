function simParams=allMutParam(p, mutx)
% Creates a matrix of which every column is a sample of mutated parameter values
% Example: Mutation mutx changes parameters p1 and p2 into the following values:
% p1=[a1, a2, a3], p2=[b1, b2, b3]
% The paramVals matrix will look like this:
% [a1 a1 a1 a2 a2 a2 a3 a3 a3;
%  b1 b2 b3 b1 b2 b3 b1 b2 b3]

% Create cell with all mutated values per parameter
n=length(mutx.params);
paramVals={};

% Loop over parameters to fill cell with all mutated values
for i=1:n
    index = mutx.params(i);
    switch mutx.direct{i}
        case '+'
            %paramVals{i}=[p(index), 2*p(index), 3*p(index), 4*p(index)];
            paramVals{i}=linspace(p(index), 10*p(index), 10);
        case '-'
            %paramVals{i}=[p(index)/4, p(index)/3, p(index)/2, p(index)];
            paramVals{i}=linspace(0.1*p(index), p(index), 10);
        case '+-'
            %paramVals{i}=[p(index)/4, p(index)/3, p(index)/2, p(index), 2*p(index), 3*p(index), 4*p(index)];
            paramVals{i}=linspace(0.1*p(index), 10*p(index), 20);
        case '0'
            paramVals{i}=[0];
        case 'specific values'       % It's a declaration for the parameter
            paramVals{i}=mutx.values;
    end
end

% Create matrix with all possibilities
amposs=1;   % Amount of possibilities, start at 1
for i=1:n
    amposs=amposs*length(paramVals{i});
end
simParams=zeros(n,amposs);

% Looping goes as follows:
% Similar to binary counting, but values up to amount of params
% Example: 2 params, with amount of values 3, 2 respectively
% The following array represents which index at the parameter we are using
% [1 1], [1 2], [2, 1], [2, 2], [3, 2]
indices=repmat(1, 1, n);
i=1;    % Simulation number

for i=1:amposs
    poss=[];
    for j=1:n
        poss=[poss, paramVals{j}(indices(j))];
    end
    simParams(:,i)=poss;

    % Increment indices: Only increase second-last param if all of the last have been looped over
    % Only increase third-last param if second-last and last param have been looped over etc.
    done=0; k=n;
    while ~done
        if indices(k)<length(paramVals{k})
            %[indices(k) length(paramVals{k})];
            indices(k)=indices(k)+1;
            done=1;
        else
            indices(k)=1;
            k=k-1;
            if k==0
                done=1;
            end
        end
    end
end
end
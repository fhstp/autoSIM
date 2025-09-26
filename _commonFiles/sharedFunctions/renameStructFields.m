function [Struct] = renameStructFields(Struct)
% This function just renames the struct fields for the variables and
% deletes the side tags in the fields. This makes life easier for plotting
% when you want to switch from the left to the right side.

% Loop fpr side
side = fieldnames(Struct);
for i_side = 1 : length(side)
    sideTmp = side{i_side};
    fnSide = fieldnames(Struct.(sideTmp));

    % Loop for Variable Category (ContacForce, ContactPressure, etc.)
    for i_vCat = 1 : length(fnSide)
        % current fieldname
        fn_vCat = fnSide{i_vCat};
        S = Struct.(sideTmp).(fn_vCat);
        if ~isempty(S)        
        % Loop per variable
        for i_Var = 1 : length(fieldnames(Struct.(sideTmp).(fn_vCat)))

            % current fieldname
            fn_vars = fieldnames(Struct.(sideTmp).(fn_vCat));
            fn_cVar = fn_vars{i_Var};

            newName = replace(fn_cVar,{['_', sideTmp], ['_', sideTmp,'_']},'');
            S = RenameField(S, fn_cVar, newName);
        end
        Struct.(sideTmp).(fn_vCat) = S;
        end
    end
end

%% Clear variables except output to prevet memory leak.
clearvars -except Struct
end

latex_dir = "../../Tesi_Andrea_Grande/Chapter_2/";
subdirs = [ ...
    "geometric_characteristics_of_a_straight_wing/", ...
   "geometric_characteristics_of_a_straight_and_tapered_wing/",...
   "geometric_characteristics_of_a_wing_with_arrow_not_null",...
    ];

% create folders (skipped if folders exist)
for i=1:length(subdirs)
    [status, msg] = mkdir(strcat(latex_dir, subdirs(i)));
end

for i=1:length(subdirs)
    status = copyfile( ...
        strcat(subdirs(i),"data.tex"), ...
        strcat(latex_dir, subdirs(i), "data.tex") ...
    );
    status = copyfile( ...
        strcat(subdirs(i),"dataset.txt"), ...
        strcat(latex_dir, subdirs(i), "dataset.txt") ...
    );
end


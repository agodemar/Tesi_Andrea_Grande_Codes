latex_dir = "../../Tesi_Andrea_Grande/Chapter_2/";
subdirs = [ ...
    "geometric_characteristics_of_a_straight_wing/", ...
   
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
end


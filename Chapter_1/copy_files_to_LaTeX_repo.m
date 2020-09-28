latex_dir = "../../Aircraft_Stability_Control_Book/";
subdir = "Chapter_1/Matlab/";
files = [ ...
    "test_0.m", "test_0.pdf", ...
    "test_1.m", "test_1.pdf", ...
    ];

for i=1:length(files)
    status = copyfile(files(i), strcat(latex_dir, subdir, files(i)));
end

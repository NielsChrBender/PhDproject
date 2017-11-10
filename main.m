%% main center %%%
close all
fprintf('------------------- Main initialized -------------------- \n');
directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m'))) % use the controlcenter.m to decide control the program
if call.fluidforce
    run(sprintf('%s%s',directory,strcat('fluidforce\main_fluidforce.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m'))) % use the controlcenter.m to decide control the program
if call.fluidforce_sphere
    run(sprintf('%s%s',directory,strcat('fluidforce\sphere\ArticlePlot_simpleSphereExample.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.contour
    run(sprintf('%s%s',directory,strcat('addons\main_contour.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.stiction
    run(sprintf('%s%s',directory,strcat('stiction\main_stiction.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.squeeze
    run(sprintf('%s%s',directory,strcat('stiction\main_squeeze.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.CFD
    run(sprintf('%s%s',directory,strcat('CFD\main_CFD.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.FEM
    run(sprintf('%s%s',directory,strcat('FEM\main_FEM.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.lifetime
    run(sprintf('%s%s',directory,strcat('lifetime\main_lifetime.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.dynamics
    run(sprintf('%s%s',directory,strcat('dynamic\main_dynamics.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.validation
    run(sprintf('%s%s',directory,strcat('validation\main_validfluidforce.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.validation_ASME17
    run(sprintf('%s%s',directory,strcat('validation\main_validfluidforce_ASME17.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.experimental_plot
    run(sprintf('%s%s',directory,strcat('validation\main_experimental.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.Leakage_study
    run(sprintf('%s%s',directory,strcat('LeakageStudy\main_leakageV3.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.single_piston
    run(sprintf('%s%s',directory,strcat('dynamic\main_single_piston.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.transient_flow || call.open_transient_flow
    run(sprintf('%s%s',directory,strcat('dynamic\main_single_piston.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.CFD_results
    figures.convergence = 0;
    figures.convergence1 = 1;
    run(sprintf('%s%s',directory,strcat('constants.m')));
    run(sprintf('%s%s',directory,strcat('dynamic\Fluidplot.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.fatigue
    run(sprintf('%s%s',directory,strcat('lifetime\main_fatigue.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.HALT
    run(sprintf('%s%s',directory,strcat('constants.m')));
    run(sprintf('%s%s',directory,strcat('lifetime\HALT.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.chamber_damping
    run(sprintf('%s%s',directory,strcat('lifetime\main_end_damping.m')))
end
%%% end of subroutine %%%

directory = strcat('C:\MAIN\'); % define your directory
run(sprintf('%s%s',directory,strcat('controlcenter.m')))
if call.neural_control
    run(sprintf('%s%s',directory,strcat('lifetime\main_end_damping.m')))
end
%%% end of subroutine %%%

fprintf('----------------------- Main end ------------------------ \n \n');
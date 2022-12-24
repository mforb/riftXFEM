try 
  testcase_ISSM
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'ISSM'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end
clear all
close all

try 
  testcase_pressure
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'pressure'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end
%try 
  %testcase_op_pen
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'pressure'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end
%clear all
%close all
%try 
  %testcase_melange_noqf_m1
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'melange 1'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end
%clear all
%close all

%try 
  %testcase_melange_noqf_m2
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'melange 2'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end
%clear all 
%close all

%%try 
  %%testcase_melangeb_m1
%%catch
  %%opfile = fopen('allrun.log','a')
  %%erm = 'Something went wrong'
  %%run = 'melange 1'
  %%fprintf(opfile,[run,'\n'])
  %%fprintf(opfile,[erm,'\n'])
  %%fclose(opfile);
%%end
%%clear all
%%close all

%try 
  %testcase_melangeb_m2
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'melange 2'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end
%clear all 
%close all

try 
  testcase_op_m1
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'melange 1 with pressure'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end
clear all
close all

%try 
  %testcase_op_m3
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'melange 2 with pressure'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end

%clear all
%close all

try 
  testcase_op_m2
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'melange 2 with pressure'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end

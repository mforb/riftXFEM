%try 
  %testcase_ISSM
%catch
  %opfile = fopen('allrun.log','a')
  %erm = 'Something went wrong'
  %run = 'ISSM'
  %fprintf(opfile,[run,'\n'])
  %fprintf(opfile,[erm,'\n'])
  %fclose(opfile);
%end

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

try 
  testcase_melange_m1
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'melange 1'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end

try 
  testcase_melange_m2
catch
  opfile = fopen('allrun.log','a')
  erm = 'Something went wrong'
  run = 'melange 2'
  fprintf(opfile,[run,'\n'])
  fprintf(opfile,[erm,'\n'])
  fclose(opfile);
end

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

%TEST_EXPM_NEW   Simple test of EXPMV_NEW.

n = 10;

for i = 1:3

   switch i
      case 1, A = gallery('invol',n);
      case 2, A = gallery('chebspec',n);
      case 3, A = [1 1e9; 0 1];
   end

   X0 = expm_new(A);
   X1 = expm_new(A,1);
   X2 = expm_new(A,2);

%    F = expm_x(A);
%    norm(F-X0,1)/norm(F,1)
%    norm(F-X1,1)/norm(F,1)
%    norm(F-X2,1)/norm(F,1)

end

disp('Test successfully completed.')

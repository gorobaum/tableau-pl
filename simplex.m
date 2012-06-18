function [ind, x] = simplex(A,b,c,m,n,print)
  [ind, x] = preparesimplex(A,b,c,m,n,print);
endfunction

function [ind, x] = preparesimplex(A,b,c,m,n,print)
  auxb = [];
  for i = 1:m
    if ( b(i) < 0.0 )
      auxb = [auxb;-b(i)];
      auxA(i,:) = -A(i,:);
    else
      auxb = [auxb;b(i)];
      auxA(i,:) = A(i,:);
    endif
  endfor
  cf1 = [];
  for i = 1:n
    cf1 = [cf1;0.0];
  endfor
  for j = (i+1):m+n
    indbase(j-i) = j;
    cf1 = [cf1;1.0];
  endfor
  
  I = eye(m);
  auxA = [auxA,I];
  auxc(1) = -(base(cf1,indbase))'*auxb;
  auxc = [auxc(1),cf1'-(base(cf1,indbase))'*auxA];
  auxA = [auxb,auxA];
  auxA = [auxc;auxA];
  
  if ( print == true )
    disp("\nSimplex: Fase 1\n")
  endif
  [ind, x, B, indb, m] = runsimplex(auxA,auxb,auxc,indbase,m,m+n,print);
  if (ind == -1 || B(1,1) != 0.0)
    ind = 1;
  else
    B = B(:,1:n+1);
    [B, indb, m] = removeslackformbase(B,m,n,indb);
    if ( print == true )
      disp("\nSimplex: Fase 2\n")
    endif
    auxc(1) = -(base(c,indb))'*B(2:m+1,1)
    auxc = [auxc(1),c'-(base(c,indb))'*B(2:m+1,2:n+1)];
    B(1,:) = auxc;
    [ind, x, B, indb] = runsimplex(B,b,c,indb,m,n,print);
  endif
  printf ("\n");
  if ( ind == -1 )
      printf ("O problema é ilimitado.\n");
    elseif (ind == 1 )
      printf ("O problema é inviável.\n");
    else
      printf ("O problema é viável e a solução ótima tem custo: %.3f\n", B(1,1));
      printf ("x = \n");
      j = 1;
      for i = 1:n
        if ( i == indb(j) )
          printf ("\t%.3f\n", B(j+1,1));
          j++;
        else
          printf ("\t0.000\n");
        endif
      endfor
    endif
endfunction

function iteration(B,indbase,j,l,iter,m,n)
  printf ("\nIteração %d\n",iter);
  printf ("\n");
  printf ("\t\t");
  for i = 1:n
    printf ("|x%d\t", i);  
  endfor
  printf ("\n");
  printf ("\t");
  aux = printf ("%.3f", B(1,1));
  if ( aux < 8 ) 
    printf ("\t");
  endif
  for i = 2:length(B(1,:))
    aux = printf ("|%.3f", B(1,i));
    if ( aux < 8 ) 
      printf ("\t");
    endif
  endfor
  printf ("\n");
  printf ("\t");
  for i = 1:8*(n+1)
    printf("-");
  endfor
  z = 1;
  for k = 2:length(B(:,1))
    printf ("\n");
    printf ("x%d\t", indbase(z));
    z++;
    printf ("%.3f", B(k,1));
    if ( aux < 8 ) 
      printf ("\t");
    endif
    for i = 2:length(B(1,:))
      if ( l != 0 && k == l && i == j )
        aux = printf ("|%.3f*", B(k,i));
        if ( aux < 8 ) 
          printf ("\t");
        endif
      else
        aux = printf ("|%.3f", B(k,i));
        if ( aux < 8 ) 
          printf ("\t");
        endif
      endif
    endfor
  endfor
  printf ("\n");
endfunction

function [ind, x, B, indb, m] = runsimplex(A,b,c,indbase,inm,n,print)
  iter = 0;
  stop = false;
  while (!stop)
    iter++;
    j = 2;
    m = inm;
    l = 0;
    while ( j <= length(A(1,:)) && A(1,j) >= 0.0 )
      j++;
    endwhile
    if ( j == length(A(1,:))+1 )
      stop = true;
      ind = 0;
    else
      u = A(:,j);
      ratio = inf;
      for i = 2:length(u)
        if ( u(i) > 0.0 && ratio > A(i,1)/u(i) )
          l = i;
          ratio = A(i,1)/u(i);
        endif
      endfor
      if ( ratio == inf)
        stop = true;
        ind = -1;
      else
        indbase(l-1) = j-1;
        if ( print == true)
          iteration(A,indbase,j,l,iter,m,n);
        endif
        for i = 1:length(u)
          if ( i != l && u(i) != 0.0 )
            A(i,:) -= (u(i)/u(l))*A(l,:);
          endif
        endfor
        A(l,:) /= u(l);
      endif
    endif
  endwhile
  indb = indbase;
  x = A(2:m+1,1);
  B = A;
  if ( print == true )
    iteration(A,indbase,j,l,iter,m,n);
  endif
endfunction

function [B, ind, m] = removeslackformbase(A,inm,n,indbase)
  B = A;
  ind = indbase;
  m = inm;
  for i = 2:m
    for j = m+1:m+n
      if (indbase(i-1) == j )
        k = 1;
        while ( k <= length(A(i,:)) && A(i,k) == 0.0 )
          k++;
        endwhile
        if ( k == length(A(i+1,:))+1 )
          B = [B(1:i-1,:);B(i+1:m+1,:)];
          ind = [ind(1:i-2),ind(i:m)];
          m--;
        else
          u = A(:,j);
          for z = 1:length(u)
            if ( z != i && u(z) != 0.0 )
              A(z,:) -= (u(z)/u(i))*A(i,:);
            endif
          endfor
          A(i,:) /= u(z);
          disp("Fodeos");
        endif
      endif
    endfor
  endfor
endfunction

function ret = base(vec, b)
  j = 1;
  ret = [];
  for i = 1:length(vec)
    if ( j <= length(b) && i == b(j) )
      ret = [ret;vec(i)];
      j++;
    endif
  endfor
endfunction

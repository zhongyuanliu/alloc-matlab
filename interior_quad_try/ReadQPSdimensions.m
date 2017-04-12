function qpDim = ReadQPSdimensions(filename)
%??????????????????????????????????????????????????????????
% ReadQPSdimensions .m
%
% This f u n c t i on s c ans a QPS f i l e f o r the s i z e s o f the QP problem and
% g i v e s
% the f o l l owi n g output .
%
% m: number o f c o n s t r a i n t s o f the form
% A( i , : ) x >= rhs ( i ) OR
% A( i , : ) x <= rhs ( i ) OR
% A( i , : ) x = rhs ( i )
%
% n : number o f v a r i a b l e s
%
% nz : the number o f nonz e ros in A
%
% qnz : the number o f o f f?diagonal e n t r i e s in the lowe r t r i a n g u l a r
% par t
% o f Q
%
% Thomas Reslow Kr¡§uth , s021898
%??????????????????????????????????????????????????????????
fid=fopen(filename, 'r' ) ;
% Read m
strLine = fgetl(fid ) ;
qpName = strLine(15 : end ) ;
strLine = fgetl(fid) ;
qpsKeyword = 0 ;
i =0;
while (~qpsKeyword)
strLine = fgetl( fid ) ;
qpsKeyword = strcmp ( strLine( 1 : 6 ) , 'COLUMN' ) ;
if (qpsKeyword==0)
i=i +1;
end
end
m = i -1;
fclose(fid) ;
fid=fopen ( filename, 'r' ) ;
strLine = fgetl( fid ) ;
qpName = strLine ( 15 : end ) ;
strLine = fgetl ( fid ) ;
qpsKeyword = 0 ;
constraintName = zeros (m+1 ,8) ;
i =0;
while (~qpsKeyword)
strLine = fgetl ( fid ) ;
qpsKeyword = strcmp ( strLine ( 1 : 6 ) , 'COLUMN' ) ;
if ( qpsKeyword==0)
i=i +1;
cN = strLine ( 5 : end ) ;
if ( length (cN)==8)
cN = cN;
elseif ( length (cN)==7)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==6)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==5)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==4)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==3)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==2)
cN = [ cN, ' ' ] ;
elseif ( length (cN)==1)
cN = [ cN, ' ' ] ;
end
end
constraintName ( i , 1 : 8 ) = cN;
end
% Read n and nz
qpsKeyword = 0 ;
n=0;
nz=0;
strold= 'C??????0';
strLine = fgetl ( fid ) ;
while (~qpsKeyword)
strnew = strLine ( 5 : 12 ) ;
if ( strcmp ( strold, strnew )==0)
n=n+1;
end
%Case : C??????1 OBJ.FUNC double R??????2 double
if ( length ( strLine )>40 && strcmp ( strLine ( 15 : 22 ) , constraintName ( 1 , 1 : 8 ))==1)
nz = nz + 1 ;
%Case : C??????1 R??????1 double R??????2 double
elseif ( length ( strLine )>40)
nz = nz + 2 ;
%Case : C??????1 OBJ.FUNC double
elseif ( length ( strLine )<39 && strcmp ( strLine ( 15 : 22 ) , constraintName( 1, 1 : 8 ) )==1)
nz = nz ;
%Case : C??????1 R??????1 double
elseif ( length ( strLine )<39)
nz = nz + 1 ;
end
strold= strLine ( 5 : 12 ) ;
strLine = fgetl ( fid ) ;
qpsKeyword = strcmp ( strLine ( 1 : 3 ) , 'RHS' ) ;
end
% Read qnz
strLine = fgetl ( fid ) ;
qpsKeyword = strcmp ( strLine ( 1 : 6 ) , 'QUADOB' ) ;
while (~qpsKeyword)
strLine = fgetl ( fid ) ;
qpsKeyword = strcmp ( strLine ( 1 : 6 ) , 'QUADOB' ) ;
end
strLine = fgetl ( fid ) ;
qpsKeyword = strcmp ( strLine ( 1 : 6 ) , 'ENDATA' ) ;
qnz=0;
while (~qpsKeyword)
strLine = fgetl ( fid ) ;
if ( length ( strLine )<39)
qnz = qnz + 1 ;
elseif ( length ( strLine )>40)
qnz = qnz + 2 ;
end
qpsKeyword = strcmp ( strLine ( 1 : 6 ) , 'ENDATA' ) ;
end
%qnz = qnz ? n ;
fclose ( fid ) ;
qpDim = struct(...
'm' , m,...
'n ' , n ,...
' nz ' , nz ,...
' qnz ' , qnz ) ;
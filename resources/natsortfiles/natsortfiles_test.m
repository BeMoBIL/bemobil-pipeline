function natsortfiles_test()
% Test function for NATSORTFILES.
%
% (c) 2014-2020 Stephen Cobeldick
%
% See also NATSORTFILES TESTFUN NATSORT_TEST NATSORTROWS_TEST
fun = @natsortfiles;
chk = testfun(fun);
%
%% Examples HTML %%
%
A = {'a2.txt', 'a10.txt', 'a1.txt'};
chk(A, fun, {'a1.txt', 'a2.txt', 'a10.txt'})
chk(A, fun, [], [3,1,2])
chk(A, fun, [], [], {{'a',2;'a',10;'a',1},{'.txt';'.txt';'.txt'}})
%
B = {'1.3.txt','1.10.txt','1.2.txt'};
chk(B, fun, {'1.2.txt', '1.3.txt', '1.10.txt'}, [3,1,2])
chk(B, '\d+\.?\d*', fun, {'1.10.txt', '1.2.txt', '1.3.txt'}, [2,3,1])
%
C = {'B.txt','10.txt','1.txt','A.txt','2.txt'};
chk(C, [],  'descend', fun, {'B.txt', 'A.txt', '10.txt', '2.txt', '1.txt'}, [])
chk(C, [], 'char<num', fun, {'A.txt', 'B.txt', '1.txt', '2.txt', '10.txt'})
%
P = 'natsortfiles_test';
S = dir(fullfile(P,'*.txt'));
chk({S.name}, fun, {'A_1.txt','A_1-new.txt','A_1_new.txt','A_2.txt','A_3.txt','A_10.txt','A_100.txt','A_200.txt'})
%
chk({'test_ccc.m'; 'test-aaa.m'; 'test.m'; 'test.bbb.m'}, fun,... D
    {'test.m'; 'test-aaa.m'; 'test.bbb.m'; 'test_ccc.m'}, [3;2;4;1])
chk({'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'}, fun,... E
    {'test.m'; 'test1.m'; 'test2.m'; 'test10.m'; 'test10-old.m'}, [3;5;1;4;2])
chk({'A2-old\test.m';'A10\test.m';'A2\test.m';'AXarchive.zip';'A1\test.m'}, fun,... F
    {'AXarchive.zip';'A1\test.m';'A2\test.m';'A2-old\test.m';'A10\test.m'}, [4;5;3;1;2])
%
G = {'1.23V.csv','-1V.csv','+1.csv','+NaNV.csv','1.200V.csv'};
chk(G, fun,...
    {'1.23V.csv', '1.200V.csv', '+1.csv', '+NaNV.csv', '-1V.csv'}, [1,5,3,4,2])
chk(G, '[-+]?(NaN|Inf|\d+\.?\d*)', fun,...
    {'-1V.csv', '+1.csv', '1.200V.csv', '1.23V.csv', '+NaNV.csv'}, [2,3,5,1,4])
%
%% Examples Mfile Help
%
chk({'A1\B', 'A+/B', 'A\B', 'A=/B', 'A/B'}, fun,... X
    {'A\B', 'A/B', 'A1\B', 'A+/B', 'A=/B'}, [3,5,1,2,4])
chk({'a2.txt', 'a10.txt', 'a1.txt'}, fun,... A
    {'a1.txt', 'a2.txt', 'a10.txt'}, [3,1,2])
chk({'test_new.m'; 'test-old.m'; 'test.m'}, fun,... B
    {'test.m'; 'test-old.m'; 'test_new.m'}, [3;2;1])
chk({'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'}, fun,... C
    {'test.m'; 'test1.m'; 'test2.m'; 'test10.m'; 'test10-old.m'}, [3;5;1;4;2])
chk({'A2-old\test.m';'A10\test.m';'A2\test.m';'A1archive.zip';'A1\test.m'}, fun,... D
    {'A1archive.zip';'A1\test.m';'A2\test.m';'A2-old\test.m';'A10\test.m'}, [4;5;3;1;2])
%
%% Orientation
%
chk({}, fun, {}, [], {}) % empty!
chk(cell(0,2,3), fun, cell(0,2,3), nan(0,2,3)) % empty!
chk({'1';'10';'20';'2'}, fun,...
    {'1';'2';'10';'20'}, [1;4;2;3])
chk({'2','10','8';'#','a',' '}, fun,...
    {'2','10','#';'8',' ','a'}, [1,3,2;5,6,4])
%
%% Index Stability
%
chk({'';'';''}, fun, {'';'';''}, [1;2;3], {{'';'';''},{'';'';''}})
%
U = {'2';'3';'2';'1';'2'};
chk(U, [], 'ascend', fun,...
    {'1';'2';'2';'2';'3'}, [4;1;3;5;2])
chk(U, [], 'descend', fun,...
    {'3';'2';'2';'2';'1'}, [2;1;3;5;4])
%
V = {'x';'z';'y';'';'z';'';'x';'y'};
chk(V, [], 'ascend', fun,...
    {'';'';'x';'x';'y';'y';'z';'z'},[4;6;1;7;3;8;2;5])
chk(V, [], 'descend', fun,...
    {'z';'z';'y';'y';'x';'x';'';''},[2;5;3;8;1;7;4;6])
%
W = {'2x';'2z';'2y';'2';'2z';'2';'2x';'2y'};
chk(W, [], 'ascend', fun,...
    {'2';'2';'2x';'2x';'2y';'2y';'2z';'2z'},[4;6;1;7;3;8;2;5])
chk(W, [], 'descend', fun,...
    {'2z';'2z';'2y';'2y';'2x';'2x';'2';'2'},[2;5;3;8;1;7;4;6])
%
%% Extension and Separator Characters
%
chk({'A.x3','A.x20','A.x','A.x1'}, fun,...
    {'A.x','A.x1','A.x3','A.x20'}, [3,4,1,2])
chk({'A=.z','A.z','A..z','A-.z','A#.z'}, fun,...
    {'A.z','A#.z','A-.z','A..z','A=.z'}, [2,5,4,3,1])
chk({'A~/B','A/B','A#/B','A=/B','A-/B'}, fun,...
    {'A/B','A#/B','A-/B','A=/B','A~/B'}, [2,3,5,4,1])
%
X = {'1.10','1.2'};
chk(X, '\d+\.?\d*', fun,...
	{'1.2','1.10'})
chk(X, '\d+\.?\d*', 'noext', fun,...
	{'1.10','1.2'})
%
Y = {'1.2','2.2','20','2','2.10','10','1','2.00','1.10'};
chk(Y, '\d+\.?\d*', fun,...
	{'1','1.2','1.10','2','2.00','2.2','2.10','10','20'},[7,1,9,4,8,2,5,6,3])
chk(Y, '\d+\.?\d*', 'noext', fun,...
	{'1','1.10','1.2','2','2.00','2.10','2.2','10','20'},[7,9,1,4,8,5,2,6,3])
%
%% Other Implementation Examples
%
% <https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/>
chk({'z1.txt','z10.txt','z100.txt','z101.txt','z102.txt','z11.txt','z12.txt','z13.txt','z14.txt','z15.txt','z16.txt','z17.txt','z18.txt','z19.txt','z2.txt','z20.txt','z3.txt','z4.txt','z5.txt','z6.txt','z7.txt','z8.txt','z9.txt'}, fun,...
    {'z1.txt','z2.txt','z3.txt','z4.txt','z5.txt','z6.txt','z7.txt','z8.txt','z9.txt','z10.txt','z11.txt','z12.txt','z13.txt','z14.txt','z15.txt','z16.txt','z17.txt','z18.txt','z19.txt','z20.txt','z100.txt','z101.txt','z102.txt'})
%
% <https://blog.jooq.org/2018/02/23/how-to-order-file-names-semantically-in-java/>
chk({'C:\temp\version-1.sql','C:\temp\version-10.1.sql','C:\temp\version-10.sql','C:\temp\version-2.sql','C:\temp\version-21.sql'}, fun,...
    {'C:\temp\version-1.sql','C:\temp\version-2.sql','C:\temp\version-10.sql','C:\temp\version-10.1.sql','C:\temp\version-21.sql'})
%
% <http://www.davekoelle.com/alphanum.html>
chk({'z1.doc','z10.doc','z100.doc','z101.doc','z102.doc','z11.doc','z12.doc','z13.doc','z14.doc','z15.doc','z16.doc','z17.doc','z18.doc','z19.doc','z2.doc','z20.doc','z3.doc','z4.doc','z5.doc','z6.doc','z7.doc','z8.doc','z9.doc'}, fun, ...
    {'z1.doc','z2.doc','z3.doc','z4.doc','z5.doc','z6.doc','z7.doc','z8.doc','z9.doc','z10.doc','z11.doc','z12.doc','z13.doc','z14.doc','z15.doc','z16.doc','z17.doc','z18.doc','z19.doc','z20.doc','z100.doc','z101.doc','z102.doc'})
%
% <https://sourcefrog.net/projects/natsort/>
chk({'rfc1.txt';'rfc2086.txt';'rfc822.txt'}, fun,...
    {'rfc1.txt';'rfc822.txt';'rfc2086.txt'})
%
% <https://www.strchr.com/natural_sorting>
chk({'picture 1.png','picture 10.png','picture 100.png','picture 11.png','picture 2.png','picture 21.png','picture 2_10.png','picture 2_9.png','picture 3.png','picture 3b.png','picture A.png'}, fun,...
    {'picture 1.png','picture 2.png','picture 2_9.png','picture 2_10.png','picture 3.png','picture 3b.png','picture 10.png','picture 11.png','picture 21.png','picture 100.png','picture A.png'})
%
chk() % display summary
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles_test
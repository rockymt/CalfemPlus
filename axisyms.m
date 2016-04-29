  function [es,et]=axisyms(er,ez,D,ed)


  colD=size(D,2);
  ee=zeros(colD,1);
  A=1/2*det([ones(3,1) er' ez']);
  z12=ez(1)-ez(2); z23=ez(2)-ez(3); z31=ez(3)-ez(1);
  r21=er(2)-er(1); r32=er(3)-er(2); r13=er(1)-er(3);
  r=sum(er)/3;
  B=1/(2*A)*[z23   0 z31   0 z12   0;
               0 r32   0 r13   0 r21;
			 r32 z23 r13 z31 r21 z12];
  B=[B(1,:);B(2,:);
     1.0/3/r 0 1.0/3/r 0 1.0/3/r 0;
	 B(3,:)];
  e=B*ed';
  if colD>4 ee([1 2 3 4])=e; else ee=e; end
  et=ee';
  es=(D*ee)';

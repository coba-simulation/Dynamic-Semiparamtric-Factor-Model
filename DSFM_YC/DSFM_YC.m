clear all
  Y=load('IT.txt');
 [m J]=size(Y);
 t0=1;%start date
 tn=m;%end date
 t=[1:1:10 15]; %maturities
 t=t./15;
 %constant over maturities (not orthogonal)
C=median(Y(:,1:end));
 for i=1:m
    Ynew(i,:)=Y(i,1:end)-C;
end
%%%%%%DSFM data preparation%%%%%%%
 TT=repmat(t',(tn - t0+1),1);
 R=rand((tn - t0+1)*J,1); 
 Da = reshape(kron((tn:-1:t0)', ones(J, 1))',J*(tn-t0+1),1);
 YY= reshape(Ynew(t0:tn,1:end)',J*(tn-t0+1),1);

 XY=[R, TT, YY, Da];
%%%%%%DSFM parameters%%%%%%%%%%%%%%

knotsX2=[t 1.01];
 knotsX1=[0 1.01];
 L=1;
 
 kX2=1;
 kX1=1;
 knotpointX1=1;
 knotpointX2=2;
epsilon = 0.0001;
maxiter = 30;
iii=10;
z0 = [rand(tn-t0+1,L+1)];

%%%%%%%DSFM estimation  without orthonormal constant factor%%%%%%%%%%%%%
  YC=NCspDSFM(XY, L, knotsX1, kX1, knotsX2, kX2, z0, epsilon, maxiter,knotpointX1,knotpointX2);
  
%%%%Goodness of fit measures%%%%%%%%%%%%  
 RMSE_IT=(mean(mean((YC.Yhat-Ynew(tn:-1:t0,:)).^2)))^.5
 EV=1-YC.description.explainedvariance(end,1)
 %%%%%%%%%%%%%%%%%%%%%%%% output plots
 figure1 = figure('PaperSize', [29.98404194812 29.67743169791], 'Visible', 'on'); 
 subplot1 = subplot(1,3,1,'Parent',figure1,'XTickLabel',{'0','5','10'},'XTick',[0 5 10],'FontSize',14,'FontWeight','bold','FontSize',16); axis on;
 subplot2 = subplot(1,3,2,'Parent',figure1,'XTickLabel',{'0','5','10'},'XTick',[0 5 10],'FontSize',14,'FontWeight','bold','FontSize',16); axis on;
 subplot3 = subplot(1,3,3,'Parent',figure1,'XTickLabel',{'0','5','10'},'XTick',[0 5 10],'FontSize',14,'FontWeight','bold','FontSize',16); axis on;

  
ylim([ceil(min(YC.mhat(1,:,1))), ceil(max(YC.mhat(1,:,1)))]);

% % Create xlabel
xlabel('Maturity','FontWeight','bold','FontSize',16);
%   title('Factor1', 'FontWeight', 'bold', 'FontSize', 16);
     plot(t.*15,YC.mhat(1,1:end-1,1),'Parent',subplot1,'Linewidth',2.5); 
     clear ylim
     ylim([ceil(min(YC.mhat(1,:,2))), ceil(max(YC.mhat(1,:,2)))]);
     xlabel('Maturity','FontWeight','bold','FontSize',16);
 plot(t.*15,YC.mhat(1,1:end-1,2),'Parent',subplot2,'Linewidth',2.5); 
  xlabel('Maturity','FontWeight','bold','FontSize',16);
         clear ylim
     ylim([ceil(min(C)), ceil(max(C))]);
  plot(t.*15,C,'Parent',subplot3,'Linewidth',2.5); 
  

  figure
  plot(YC.z,'Linewidth',2.5)
  xlabel('Time','FontWeight','bold','FontSize',16);
 
 
 
%  %%%%%%%%factors
%  plot(t.*15,YC.mhat(1,1:end-1,2),'Linewidth',2.5); 
%  hold on
%  plot(t.*15,YC.mhat(1,1:end-1,1),'Linewidth',2.5,'Color','red'); 
%   hold on
%  plot(t.*15,C,'Linewidth',2.5,'Color','green'); 
%  %%%%%%%factor loadings
%  figure 
%  plot(YC.z(1:end,1:2))
 
 

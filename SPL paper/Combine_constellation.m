clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The QPSK modulation %%%  

tic;        %Timer Started

%Two modulation Schemes QPSK & QAM
% QPSK <--- Z1
% QAM16<--- Z2

l= floor(100000/12)*12;            %No. of Bits @Note: should be divisble by nbits

snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/20);
berc1 = zeros(1,length(snrdb));
berc2 = zeros(1,length(snrdb));


bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zpsk16 = [1 cos(pi/8)+sin(pi/8)*1i cos(3*pi/8)+sin(3*pi/8)*1i cos(pi/4)+sin(pi/4)*1i -cos(pi/4)+sin(pi/4)*1i -cos(pi/8)+sin(pi/8)*1i -cos(3*pi/8)+sin(3*pi/8)*1i 1i -1i cos(3*pi/8)-sin(3*pi/8)*1i cos(pi/8)-sin(pi/8)*1i cos(pi/4)-sin(pi/4)*1i -cos(pi/4)-sin(pi/4)*1i -cos(3*pi/8)-sin(3*pi/8)*1i -cos(pi/8)-sin(pi/8)*1i -1];    %psk16 Grey coding "Unit Power"
zqpsk = (1/sqrt(2))*[1+1i -1+1i 1-1i -1-1i];    %qpsk Grey coding "Unit Power"

M1 = length(zqpsk);
M2 = length(zpsk16);

nbits1 = log2(M1);
sym1=zeros(1,l/nbits1);   %Symbol sent -> QPSK

nbits2 = log2(M2);
sym2=zeros(1,l/(nbits2+nbits1));   %Symbol sent -> Qpsk16[ -|-| -|-|-|- ]

for len = nbits1:nbits1:l
    switch (bi2de([bits(len) bits(len-1)]))
        case 0
            sym1(len/nbits1) = zqpsk(1);
        case 1
            sym1(len/nbits1) = zqpsk(2);         
        case 2
            sym1(len/nbits1) = zqpsk(3);            
        case 3
             sym1(len/nbits1) = zqpsk(4);           
    end
end


for len = nbits2+2:nbits2+2:l
    switch (bi2de([bits(len) bits(len-1) bits(len-2) bits(len-3)]))
        case 0
            sym2(len/(nbits2+nbits1))= zpsk16(1);
        case 1
            sym2(len/(nbits2+nbits1))= zpsk16(2);         
        case 2
            sym2(len/(nbits2+nbits1))= zpsk16(3);            
        case 3
            sym2(len/(nbits2+nbits1))= zpsk16(4);           
        case 4
            sym2(len/(nbits2+nbits1))= zpsk16(5);            
        case 5
            sym2(len/(nbits2+nbits1))= zpsk16(6);            
        case 6
            sym2(len/(nbits2+nbits1))= zpsk16(7);            
        case 7
            sym2(len/(nbits2+nbits1))= zpsk16(8);            
        case 8
            sym2(len/(nbits2+nbits1))= zpsk16(9);        
        case 9
            sym2(len/(nbits2+nbits1))= zpsk16(10);            
        case 10
            sym2(len/(nbits2+nbits1))= zpsk16(11);            
        case 11
            sym2(len/(nbits2+nbits1))= zpsk16(12);          
        case 12
            sym2(len/(nbits2+nbits1))= zpsk16(13);            
        case 13
            sym2(len/(nbits2+nbits1))= zpsk16(14);            
        case 14
            sym2(len/(nbits2+nbits1))= zpsk16(15);            
        case 15            
            sym2(len/(nbits2+nbits1))= zpsk16(16);           
    end
end

out = zeros(1,l/nbits1);     %Received symbols array
outc1 = zeros(1,(2*l)/(nbits1+nbits2));

SD = 1/sqrt(2);
loop_Cond = sym1;

for Case=1:2
    for n=1:length(snrdb)
        if Case==1
            loop_Cond = size(sym1,2);
        elseif Case==2
            loop_Cond = size(sym2,2)*2;
        end
        
        for loop=1:2:loop_Cond        %Two symbols at every iteration.

            Hx11=normrnd(0,SD,[1 1]);
            Hy11=normrnd(0,SD,[1 1]);
            H1=(Hx11+1i*Hy11);

            Hx22=normrnd(0,SD,[1 1]);
            Hy22=normrnd(0,SD,[1 1]);
            H2=(Hx22+1i*Hy22);

            w1_0 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RH
            w1_1 = snrx(n)*(1/sqrt(2))*(randn(1,1) + 1i*randn(1,1));        %Noise Profile RV

            if Case==1
                q1=sym1(loop);                   % First symbol
                q2=sym1(loop+1);                 % Second symbol
            elseif Case==2
                temp_loop = floor(loop/2)+1;
                q1=sym1(temp_loop+(temp_loop-1)*2);                   % First symbol
                
                q2=sym2(temp_loop);                              % Second symbol
            end
            


            RH = q1*H1 - q2*conj(H2) + w1_0;      %Received at RH            
            RV = q1*H2 + q2*conj(H1) + w1_1;      %Received at RV

    %%%%%%%% Received Signal at Decoder %%%%%%%%%%%

            if Case==1
                lenz1 = size(zqpsk);
                lenz2 = size(zqpsk);
            elseif Case==2
                lenz1 = size(zqpsk);
                lenz2 = size(zpsk16);
            end
            
            const=-100000;
            for z1=1:lenz1(2)
                for z2=1:lenz2(2)

                    C1 = zqpsk(z1);
                    
                    if Case==1  
                        C2 = zqpsk(z2);
                    elseif Case==2
                        C2 = zpsk16(z2);
                    end

                    DR = (real(conj(RH)*(C1*H1-C2*conj(H2)) + RV*(conj(C1)*conj(H2)+conj(C2)*H1)));       %Decoding Symbol
                    if (DR > const)
                        temp = [z1 z2];
                        const= DR;
                    end
                end
            end
            
            

            if Case==1  
                out(loop) = zqpsk(temp(1));
                out(loop+1) = zqpsk(temp(2));
            elseif Case==2
                outc1(loop) = zqpsk(temp(1));
                outc1(loop+1) = zpsk16(temp(2));
            end            
        end
        
        if Case==1  
            for k=1:1:l/nbits1
                switch out(k)
                    case zqpsk(1)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 0];
                    case zqpsk(2)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[0 1];                
                    case zqpsk(3)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 0];                
                    case zqpsk(4)
                       bit(nbits1*k-(nbits1-1):nbits1*k)=[1 1];
                end
            end
            berc1(n)=(l-sum(bits==bit))/l;
        elseif Case==2
            for k=1:2:2*l/(nbits1+nbits2)
                ktemp = floor(k/2)+1;
                switch outc1(k)
                    case zqpsk(1)
                       bit((nbits1+nbits2)*ktemp-((nbits1+nbits2)-1):((nbits1+nbits2)*ktemp)-nbits2)=[0 0];
                    case zqpsk(2)
                       bit((nbits1+nbits2)*ktemp-((nbits1+nbits2)-1):((nbits1+nbits2)*ktemp)-nbits2)=[0 1];                
                    case zqpsk(3)
                       bit((nbits1+nbits2)*ktemp-((nbits1+nbits2)-1):((nbits1+nbits2)*ktemp)-nbits2)=[1 0];                
                    case zqpsk(4)
                       bit((nbits1+nbits2)*ktemp-((nbits1+nbits2)-1):((nbits1+nbits2)*ktemp)-nbits2)=[1 1];
                end
                
                switch outc1(k+1)
                    case zpsk16(1)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 0 0 0];
                    case zpsk16(2)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 0 0 1];                
                    case zpsk16(3)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 0 1 0];                
                    case zpsk16(4)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 0 1 1];
                    case zpsk16(5)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 1 0 0];               
                    case zpsk16(6)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 1 0 1];                
                    case zpsk16(7)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 1 1 0];                
                    case zpsk16(8)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[0 1 1 1];                
                    case zpsk16(9)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 0 0 0];               
                    case zpsk16(10)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 0 0 1];                
                    case zpsk16(11)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 0 1 0];                
                    case zpsk16(12)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 0 1 1];                
                    case zpsk16(13)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 1 0 0];               
                    case zpsk16(14)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 1 0 1];                
                    case zpsk16(15)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 1 1 0];                
                    case zpsk16(16)
                       bit((nbits1+nbits2)*ktemp-(nbits2-1):(nbits1+nbits2)*ktemp)=[1 1 1 1];        
                end
            end
            berc2(n)=(l-sum(bits==bit))/l;
        end  
    end
end
    
semilogy(snrdb,berc1,'m-*','linewidth',2)
hold on
semilogy(snrdb,berc2,'r-*','linewidth',2)
grid on
legend('QPSK,QPSK','QPSK,PSk16')
xlim([-4 30])
ylim([10^-8 10^0])
title('1x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 


        
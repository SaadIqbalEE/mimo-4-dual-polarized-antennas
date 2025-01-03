clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The QPSK modulation %%%  

tic;        %Timer Started

%Two modulation Schemes QPSK & QPSKr
% QPSK <--- Z1
% QPSKr<--- Z2

l= floor(100000/8)*8;            %No. of Bits @Note: should be divisble by nbits

snrdb=-4:2:30;      %Array of SNR
snrx= 10.^(-snrdb/20);
ber4c1 = zeros(1,length(snrdb));
ber4c2 = zeros(1,length(snrdb));


bits=round(rand(1,l));  %Bits generation
bit = 0;                %Resultant bits

zqpsk = exp((1i*(pi/2)*[0:3])+(0));          %qpsk "Unit Power"
zqpskr  = exp((1i*(pi/2)*[0:3])+(pi/4));       %qpsk rotated "Unit Power"

M1 = length(zqpsk);
M2 = length(zqpskr);

nbits1 = log2(M1);
sym1=zeros(1,l/nbits1);   %Symbol sent -> QPSK

nbits2 = log2(M2);
sym2=zeros(1,l/nbits2);   %Symbol sent -> QPSK

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

for len = nbits2:nbits2:l
    switch (bi2de([bits(len) bits(len-1)]))
        case 0
            sym2(len/nbits1) = zqpskr(1);
        case 1
            sym2(len/nbits1) = zqpskr(2);         
        case 2
            sym2(len/nbits1) = zqpskr(3);            
        case 3
            sym2(len/nbits1) = zqpskr(4);           
    end
end

out = zeros(1,l/nbits1);     %Received symbols array

SD = 1/sqrt(2);
loop_Cond = sym1;

for Case=1:2
    for n=1:length(snrdb)
        for loop=1:2:size(sym1,2)        %Two symbols at every iteration.

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
                q1=sym1(loop);                   % First symbol
                q2=sym2(loop+1);                 % Second symbol
            end

            RH = q1*H1 - q2*conj(H2) + w1_0;      %Received at RH            
            RV = q1*H2 + q2*conj(H1) + w1_1;      %Received at RV

    %%%%%%%% Received Signal at Decoder %%%%%%%%%%%

            if Case==1
                lenz1 = size(zqpsk);
                lenz2 = size(zqpsk);
            elseif Case==2
                lenz1 = size(zqpsk);
                lenz2 = size(zqpskr);
            end
            
            const=-100000;
            for z1=1:lenz1(2)
                for z2=1:lenz2(2)

                    C1 = zqpsk(z1);
                    
                    if Case==1  
                        C2 = zqpsk(z2);
                    elseif Case==2
                        C2 = zqpskr(z2);
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
                out(loop) = zqpsk(temp(1));
                out(loop+1) = zqpskr(temp(2));
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
            ber4c1(n)=(l-sum(bits==bit))/l;
        elseif Case==2
            for k=1:1:l/nbits1
                if (rem(k,2) == 1)
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
                else
                    switch out(k)
                        case zqpskr(1)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 0];
                        case zqpskr(2)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[0 1];                
                        case zqpskr(3)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 0];                
                        case zqpskr(4)
                           bit(nbits2*k-(nbits2-1):nbits2*k)=[1 1];
                    end
                end
            end
            ber4c2(n)=(l-sum(bits==bit))/l;
        end  
    end
end
    
semilogy(snrdb,ber4c1,'m-*','linewidth',2)
hold on
semilogy(snrdb,ber4c2,'r-*','linewidth',2)
grid on
legend('QPSK,QPSK','QPSK,QPSKr')
xlim([-4 30])
ylim([10^-8 10^0])
title('1x1 system BER vs SnR Curve');
xlabel(' SNR (dB)') % x-axis label
ylabel(' BER ') % y-axis label
toc 


        
%questions collection

m = 8;
N = 255;
K = 239;
numPacket=9;

puncturepattern12=[1;1];                                 %puncture pattern when the code rate=1/2
puncturepattern23=[1;1;0;1];                                                                                                   %2/3
puncturepattern56=[1;1;0;1;1;0;0;1;1;0];                                                                                %5/6
puncturepatternSelect=[puncturepattern12;puncturepattern23;puncturepattern56]; % a column array (2+4+10)*1 
PC=[1,2;3,6;7,16];           % work as a point to indicate the location(start,end position) of different puncture pattern. 
                                            %e.g. The puncture pattern for 2/3 start at the No 3 element, end at No 6
                                            %element of the array'puncturepatternSelect', so the pattern is [1;0;1;1]
tracebackdepthSelect=[60,75,120];  %3 tracebackdepths for 3 code rates
coderateSelect=[1/2,2/3,5/6];
M=4;
EbNo=-4:0.2:30;
countern=[65280 48960 39168];

%%%%%%%%%%%%%%%uncoded
EbNoun=-4:0.5:16;



for b=1:3
    

%%%%    
  EbNoratesemilogy=zeros(1,length(EbNoun)); 
 hErrorCalc = comm.ErrorRate(); 
%hError = comm.ErrorRate();
BER8=zeros(1,length(EbNoun));
 hMod = comm.PSKModulator(M, 'BitInput',true); 
 hDemod = comm.PSKDemodulator(M, 'BitOutput',true);
%hQPSKMod=comm.QPSKModulator;
%hQPSKDemod=comm.QPSKDemodulator;

hChan= comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)');
hChan.BitsPerSymbol=2;

BERVec1=zeros(3,length(EbNoun));
EbNoEncoderOutput = EbNoun +10*log10(coderateSelect(b)*188/204);
for n=1:length(EbNoEncoderOutput)
     reset(hErrorCalc);
     [EbNoEncoderOutput(n),coderateSelect(b)]
    hChan.EbNo = EbNoEncoderOutput(n);
    
     %hChanQPSK.EbNo =SNR(n);
     
     while (BERVec1(2,n)<100)&&(BERVec1(3,n)<5000000)
         data0 = randi([0,1],64800, 1);
         
         modQ=step(hMod,data0);
         chanQ=step(hChan,modQ);
         demodQ=step(hDemod,chanQ);

         BERVec1(:,n)= step(hErrorCalc, data0, demodQ);
         BER8(n)=BERVec1(1,n);
     end
end
semilogy(EbNoun, BER8,'-');
%axis([-6 12 8e-6 10]);
hold on;
grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%coded
ShortenedLength = 51;
S = K - ShortenedLength; %188
SS = N - ShortenedLength; %204
RS_Rate = S/(N - ShortenedLength);

gp = rsgenpoly(N,K,[],0); %generate polynomial
BERall=zeros(3,length(EbNo));

for gs=1:3
%gs=3;

%rayleigh channel
chanx=normrnd(0,sqrt(1/2),countern(gs),1);
chany=normrnd(0,sqrt(1/2),countern(gs),1);
hChanray=abs(chanx+chany*sqrt(-1));

BERp=zeros(3,length(EbNo));
for b=1:length(EbNo)
   %1st for loop here for different code rate!!!!!!!!

%conv code encode
hConvEnc = comm.ConvolutionalEncoder(poly2trellis(7,[133 171]));       %maybe need to change the 133 171
hConvEnc.PuncturePatternSource = 'Property';
hConvEnc.PuncturePattern =puncturepatternSelect(PC(gs,1):PC(gs,2),1);

%QPSK modulation & demodulation
hMod = comm.QPSKModulator;
hDemod = comm.QPSKDemodulator;

%AWGNchannel
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)');           %signalpower not sure

%conv code decode
hVitDec = comm.ViterbiDecoder(poly2trellis(7,[133 171]),...
    'InputFormat','Hard');
hVitDec.PuncturePatternSource  = 'Property';                           %property not sure
hVitDec.PuncturePattern = hConvEnc.PuncturePattern;
hVitDec.TracebackDepth =tracebackdepthSelect(gs) ;

% add fuction Eb/N0
EbNorate= EbNo(b)+10*log10(coderateSelect(gs));

hChan.EbNo =EbNorate;
hChan.BitsPerSymbol=2;
%ErrorCalculator
hErrorCalc = comm.ErrorRate('ReceiveDelay',hVitDec.TracebackDepth);
hError=comm.ErrorRate;

%hChan.EbNo = 30; % 2nd for loop here for different SNR!!!!!!!!

%RS encoder & decoder
enc = comm.RSEncoder(N,K,gp,S);
dec = comm.RSDecoder(N,K,gp,S);

%ConvolutionalInterleaver & ConvolutionalDeinterleaver
nrows = 12; slope = 17; % Interleaver parameters
%D = nrows*(nrows-1)*slope; % Delay of interleaver/deinterleaver pair
hInt = comm.ConvolutionalInterleaver('NumRegisters',nrows,'RegisterLengthStep', slope);
hDeint = comm.ConvolutionalDeinterleaver('NumRegisters',nrows,'RegisterLengthStep', slope);


% Calculate BER for length(SNR) frames

% while loop here for collect 100 bits!!!!!!!!

while ((BERp(2,b)<100)&&(BERp(3,b)<20000000))
 %[b BERp(2,b) gs]   
data = randi([0 255],S*numPacket,1);  %£¨188*9£©*1 bytes
data0 = zeros(S*11,1);       %2068*1 bytes
dataRSb = [data;data0];    %3760*1 bytes  1692/3760

encdata=step(enc,dataRSb);            %4080*1 bytes   RS Encode 
convintldata=step(hInt,encdata);    %4080*1 bytes   convolutional interleaver 

datainB=de2bi(convintldata)';  % 4080*8 -> 8*4080 bits  10jinzhi turn to 2jinzhi
datainBb=datainB(:);  % 32640*1 bits 

encData = step(hConvEnc, datainBb);  %convolutional code
%encData = step(hConvEnc, encdata); 
modSignal = step(hMod,encData); %QPSK modulation
raydata=hChanray.*modSignal;% Rayleigh channel
channelOutput = step(hChan,raydata); %AWGN channel
DemodData = step(hDemod, channelOutput); % QPSK demodulation
decData = step(hVitDec, (DemodData));% 32640*1 bits viterbi decoder

decDataTB = [decData((hVitDec.TracebackDepth+1):end,1); zeros(hVitDec.TracebackDepth,1)];% (32640-60+60)*1 bits
dataoutB=reshape(decDataTB,8,SS*(11+numPacket))'; %4080*8 bits
dataoutD=bi2de(dataoutB); % 4080*1 bytes 
%dataoutDm=reshape(dataoutD, 12,204)';
deconvintldata=step(hDeint,dataoutD); % 4080*1 bytes  Deinterleaver

decdata=step(dec,deconvintldata); %3760 bytes RS Decoder
deout=decdata(((nrows-1)*S+1):end,1); %1692 bytes

error1=step(hErrorCalc,datainBb,decData);
error2=step(hError,data,deout);

%BER=step(hErrorCalc, data, decData);

BERp(:,b) = error2;
%[error1(1) error2(1)]
end

 [EbNo(b) BERp(2,b) gs BERp(1,b)]
 
 if (BERp(1,b)<=2e-5)
     break
 end
end
BERall(gs,:)=BERp(1,:);
semilogy(EbNo,BERall(gs,:),'-');
hold on;

end
%legend('uncoded1/2','uncoded2/3','uncoded5/6','RS1/2','RS2/3','RS5/6');
axis([-4 28 8e-6 10]);
grid on;

legend('uncoded 1/2','uncoded 2/3','uncoded 5/6','RS1/2','RS2/3','RS5/6');
title('RS(188,204)');
xlabel('EB/N0(dB)');
ylabel('BER')
//Integration methode des trapeze simplifiee
function z_out = trapez(y)
        z_out = sum(y(1:$-1) + y(2:$));
endfunction

function fmax = TDO_matlab(inputfile, outputfile)
    // Parametres de la manip
    fourier_le = 1024;
    time_le    = 1024;
    dfmin      = 1;
    dt         = 2e-8;
    
    //stacksize(1.0D+8);
    inpf = fileinfo(inputfile)
    disp(inpf)
    datatoread = inpf(1) - 512;
    nbblocks = int((datatoread / (65536*2)) - 1) 
// Chargement des donnees
    //fid = mopen(inputfile, 'rb');
    //i = 0;
    //mseek(512,fid,'set');
    //mclearerr(fid)
    //V=[];
    //i = 1
    //while i<= nbblocks do    
    //     W=mget(65536,'s',fid);
    //     V(i,:)=W;
    //     i = i + 1;
    //end
    //V=matrix(V',int(65536*nbblocks),1)
    //mclose(fid);
    //write("v.dat",V);
    disp('start')
    fd1=mopen(inputfile,'rb');
    disp('read')
    V = mget((inpf(1)-512)/2,'s',fd1)
    V=V'
    disp(size(V))
    disp('done')
    mclose(fd1)
    
    // Definition du nombre de pas en fonction des params
    pstart     = 2;
    pend       = floor(size(V,1)/time_le)-200;
    // Soustraction frequence nulle
    V = V - mean(V);
    
    // Vecteur temps constant
    [l,c]=size(V);
    t = 0:dt:l* dt;
    
    // FFT et obtention de la frequence approximative
    tf         = 0:1/((time_le-1)*dt):1/dt;
    Vf         = abs(real(fft(V(1:time_le))));
    fmax       = zeros(pend, 2);
    [ans,pos]    = max(Vf(1:floor(time_le/2)));
    fmax(1, 2) = tf(pos);
    fmax(1, 1) = 0;

    // Calcul des constantes de la boucle
    t_le    =  t(1:fourier_le); 
    expon   = -2 * %i * %pi * t_le;
    deltaf0 =  (tf(2)-tf(1))/1000;
    
    // Fenetrage Fourier (comment useless lines)
    // Rectangulaire
    //win = ones([1,fourier_le]);
    // Cosinus
    win = 0:fourier_le-1;
    win = 1-cos(win*2*%pi/(fourier_le-1));
    // Transposee de la matrice
    win = win';
    
    disp('end loading')
    for i=pstart:pend
        // Utilisation de la dernière valeur comme point de depart
        a = fmax(i-1, 2);
        
        // On calcule la valeur par dichotomie    
        deltaf = deltaf0;
        Vtemp = V((i-1)*time_le+1:(i-1)*time_le+fourier_le).*win;
        
        // Poids spectral frequence de reference
        F = abs(trapez(Vtemp.*exp(expon'*a)));
        
        // Si poids plus grand pour frequence plus grande      
        if abs(trapez(Vtemp.*exp(expon'*(a+deltaf)))) > F
            essaimax = F;
            
            while abs(deltaf)>dfmin
                F = abs(trapez(Vtemp.*exp(expon'*(a+deltaf))));
                if F > essaimax
                    essaimax = F;
                    a = a+deltaf;
                else
                    deltaf = -deltaf/5;
                end
            end
        // Si poids plus faible pour frequence plus grande
        else
            essaimax = F;
            
            while abs(deltaf)>dfmin
                F = abs(trapez(Vtemp.*exp(expon'*(a-deltaf))));
                if F > essaimax
                    essaimax = F;
                    a = a-deltaf;
                else
                    deltaf = -deltaf/5;
                end
            end
        end
     
        // Storage
        fmax(i, 2) = a;
        fmax(i, 1) = ((i-1)*time_le+fourier_le/2)*dt;
    end
    //w = mean(fmax(2:50,2));
    w = 0.0
    // Plot results
    plot(fmax(:,1),-(fmax(:,2)-w));
    
endfunction

A = TDO_matlab('input.dat', 'output.dat');



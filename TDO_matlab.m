function fmax = TDO_matlab(inputfile, outputfile)
    % Integration methode des trapeze simplifiee
    function z_out = trapez(y)
        z_out = sum(y(1:end-1) + y(2:end));
    end

    % Parametres de la manip
    fourier_le = 1024;
    time_le    = 1024;
    dfmin      = 1;
    dt         = 2e-8;
    
    % Chargement des donnees
    fid = fopen(inputfile, 'r');
    fseek(fid,512,'bof');
    V = fread(fid, 'int16');
    fclose(fid);
    
    % Definition du nombre de pas en fonction des params
    pstart     = 2;
    pend       = floor(size(V,1)/time_le)-200;
    
    % Soustraction frequence nulle
    V = V-mean(V);
    % Vecteur temps constant
    t = 0:dt:size(V)*dt;
    
    % FFT et obtention de la frequence approximative
    tf         = 0:1/((time_le-1)*dt):1/dt;
    Vf         = abs(real(fft(V(1:time_le))));
    fmax       = zeros(pend, 2);
    [~,pos]    = max(Vf(1:floor(time_le/2)));
    fmax(1, 2) = tf(pos);
    fmax(1, 1) = 0;

    % Calcul des constantes de la boucle
    t_le    =  t(1:fourier_le); 
    expon   = -2j*pi*t_le;
    deltaf0 =  (tf(2)-tf(1))/1000;
    
    % Fenetrage Fourier (comment useless lines)
    % Rectangulaire
    window = ones([1,fourier_le]);
    % Cosinus
    window = 0:fourier_le-1;
    window = 1-cos(window*2*pi/(fourier_le-1));
    % Transposee de la matrice
    window = window';
    
    for i=pstart:pend
        % Utilisation de la dernière valeur comme point de depart
        a = fmax(i-1, 2);
        
        % On calcule la valeur par dichotomie    
        deltaf = deltaf0;
        Vtemp = V((i-1)*time_le+1:(i-1)*time_le+fourier_le).*window;
        
        % Poids spectral frequence de reference
        F = abs(trapez(Vtemp.*exp(expon'*a)));
        
        % Si poids plus grand pour frequence plus grande      
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
        % Si poids plus faible pour frequence plus grande
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
     
        % Storage
        fmax(i, 2) = a;
        fmax(i, 1) = ((i-1)*time_le+fourier_le/2)*dt;
    end
    
    % Plot results
    plot(fmax(:,1),-(fmax(:,2)-mean(fmax(2:50,2))));
    
    dlmwrite(outputfile, fmax, ' ');
end
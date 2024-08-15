function generateplotspectrum()
    val = -10:0.1:10;
    letters = ["a","b","c","d","e","f"];
    for j = 1:length(letters)
        letter = letters(j);
        switch letter
            case 'a'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[val(i) 1 1 1 1 1],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/a/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end
            case 'b'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[1 val(i) 1 1 1 1],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/b/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end
            case 'c'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[1 1 val(i) 1 1 1],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/c/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end
            case 'd'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[1 1 1 val(i) 1 1],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/d/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end


            case 'e'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[1 1 1 1 val(i) 1],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/e/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end
            case 'f'
                for i = 1:length(val)
                    CEF_MagnetotropicFittingAndCalculation('KYbSe2 18T','Display Fit Results k',[[1 1 1 1 1 val(i)],[0.803033,1.293324,0]])
                    
                    saveas(gcf,['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/plots/CEF Fits/f/',[convertStringsToChars(letter),num2str(i),'.jpg']])
                    close(gcf);
                end
        end
    end
%Plot Fixer
%Written by: Matt Svrcek  12/05/2001 Edited by: William Greenbaum 10/22/2013

%Run this script after generating the raw plots.  It will find
%all open figures and adjust line sizes and text properties.

%Change the following values to suit your preferences.  The variable
%names and comments that follow explain what each does and their options.

Plot.plotlsize = 2; %thickness of plotted lines, in points
Plot.axislsize = 1; %thickness of tick marks and borders, in points
Plot.markersize = 10;  %size of line markers, default is 6

%font names below must exactly match your system's font names
%check the list in the figure pull down menu under Tools->Text Properties
%note, the script editor does not have all the fonts, so use the figure menu

Plot.axisfont = 'Helvetica'; %changes appearance of axis numbers
Plot.axisfontsize = 16;            %in points
Plot.axisfontweight = 'normal';    %options are 'light' 'normal' 'demi' 'bold'
Plot.axisfontitalics = 'normal';   %options are 'normal' 'italic' 'oblique'

Plot.legendfont = 'Helvetica'; %changes text in the legend
Plot.legendfontsize = 15;
Plot.legendfontweight = 'normal';
Plot.legendfontitalics = 'normal';

Plot.labelfont = 'Helvetica';  %changes x, y, and z axis labels
Plot.labelfontsize = 18;
Plot.labelfontweight = 'bold';
Plot.labelfontitalics = 'normal';

Plot.titlefont = 'Helvetica';  %changes title
Plot.titlefontsize = 14;
Plot.titlefontweight = 'normal';
Plot.titlefontitalics = 'normal';

Plot.textfont = 'Helvetica';   %changes text
Plot.textfontsize = 15;
Plot.textfontweight = 'normal';
Plot.textfontitalics = 'normal';


%stop changing things below this line
%----------------------------------------------------
Plot.figureh = findobj('Type', 'figure');
Plot.axesh = findobj('Type', 'axes'); 
Plot.lineh = findobj(Plot.axesh, 'Type', 'line');
Plot.axestexth = findobj(Plot.axesh, 'Type', 'text');
set(Plot.lineh, 'LineWidth', Plot.plotlsize)
set(Plot.lineh, 'MarkerSize', Plot.markersize)
set(Plot.axesh, 'LineWidth', Plot.axislsize)
set(Plot.axesh, 'FontName', Plot.axisfont)
set(Plot.axesh, 'FontSize', Plot.axisfontsize)
set(Plot.axesh, 'FontWeight', Plot.axisfontweight)
set(Plot.axesh, 'FontAngle', Plot.axisfontitalics)
set(Plot.axestexth, 'FontName', Plot.textfont)
set(Plot.axestexth, 'FontSize', Plot.textfontsize)
set(Plot.axestexth, 'FontWeight', Plot.textfontweight)

% Edit Labels
for(i = 1:1:size(Plot.axesh))

    Plot.curAx = Plot.axesh(i);
%     (strcmpi(get(axesh(i), 'tag'), 'legend') == false);
    if (strcmpi(get(Plot.axesh(i), 'tag'), 'legend') == false)
        set(get(Plot.curAx,'XLabel'), 'FontName', Plot.labelfont)
        set(get(Plot.curAx,'XLabel'), 'FontSize', Plot.labelfontsize)
        set(get(Plot.curAx,'XLabel'), 'FontWeight', Plot.labelfontweight)
        set(get(Plot.curAx,'XLabel'), 'FontAngle', Plot.labelfontitalics)
        set(get(Plot.curAx,'YLabel'), 'FontName', Plot.labelfont)
        set(get(Plot.curAx,'YLabel'), 'FontSize', Plot.labelfontsize)
        set(get(Plot.curAx,'YLabel'), 'FontWeight', Plot.labelfontweight)
        set(get(Plot.curAx,'YLabel'), 'FontAngle', Plot.labelfontitalics)
        set(get(Plot.curAx,'ZLabel'), 'FontName', Plot.labelfont)
        set(get(Plot.curAx,'ZLabel'), 'FontSize', Plot.labelfontsize)
        set(get(Plot.curAx,'ZLabel'), 'FontWeight', Plot.labelfontweight)
        set(get(Plot.curAx,'ZLabel'), 'FontAngle', Plot.labelfontitalics)
        set(get(Plot.curAx,'Title'), 'FontName', Plot.titlefont)
        set(get(Plot.curAx,'Title'), 'FontSize', Plot.titlefontsize)
        set(get(Plot.curAx,'Title'), 'FontWeight', Plot.titlefontweight)
        set(get(Plot.curAx,'Title'), 'FontAngle', Plot.titlefontitalics)
    else
        set(Plot.curAx, 'FontName', Plot.labelfont)
        set(Plot.curAx, 'FontSize', Plot.labelfontsize)
        set(Plot.curAx, 'FontWeight', Plot.labelfontweight)
        set(Plot.curAx, 'FontAngle', Plot.labelfontitalics)
    end
    set(get(Plot.curAx,'legend'),'FontSize', Plot.legendfontsize);
end

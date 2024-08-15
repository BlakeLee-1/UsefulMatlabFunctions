intervals = length(datacell.sectionvec);

for i = 1:intervals
    if i ~= length(datacell.sectionvec)
        sections(i).Time = datacell.Time(datacell.sectionvec(i):datacell.sectionvec(i+1));
        sections(i).Voltage = datacell.Thermometer_Voltage(datacell.datacell.sectionvec(i+1));
        minn = min(sections(i).Time);
        for j = 1:length(sections(i).Time)
            sections(i).Time(j) = sections(i).Time(j) - minn;
        end
    elseif i == length(datacell.sectionvec)
        sections(i).Time = datacell.Time(datacell.sectionvec(i):end);
        sections(i).Voltage = datacell.Thermometer_Voltage(datacell.sectionvec(i):end);
        minn = min(sections(i).Time);
        for j = 1:length(sections(i).Time)
            sections(i).Time(j) = sections(i).Time(j) - minn;
        end
    end
end
colorz = lines(intervals);
figure; xlabel('Time (s)'); ylabel('Voltage (\muV)'); title('Square Wave Response'); box on; grid on;

for i = 1:intervals
    
    hold on;
    plot(sections(i).Time,sections(i).Voltage, 'o', 'Color', colorz(i,:), 'MarkerFaceColor', colorz(i,:), 'DisplayName', ['Section ', num2str(i)])
    sections(i).Time
    
end


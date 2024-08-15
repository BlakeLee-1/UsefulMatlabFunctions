function [x,y,ExtraData] = ExtractPlotData_Single(AxesHandle,Child)

x = AxesHandle.Children(Child).XData;
y = AxesHandle.Children(Child).YData;
ExtraData.Color = AxesHandle.Children(Child).Color;
ExtraData.DisplayName = AxesHandle.Children(Child).DisplayName;

end
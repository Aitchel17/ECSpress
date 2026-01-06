function analyze_shaded_graph(x_array,mean,ci,color)
fill([x_array fliplr(x_array)], [(mean+ci)' fliplr((mean-ci)')], color, 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(x_array, mean, 'color',color, 'LineWidth', 2)
end


package org.example;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.util.Arrays;



public class LeastSquares {
    public static double a0, a1, a2;

    public static void main(String[] args) {
        double[][] data = {
                {0, 0, 0, 1, 1, 2, 2, 2},
                {1.5, 2.5, 3.5, 1.5, 3.5, 1.5, 2.5, 2.5},
                {2.3, 4.9, 1.7, 4.4, 3.4, 6.7, 6.2, 7.2}
        };

        int n = data[0].length;
        double sumX1 = 0, sumX2 = 0, sumY = 0;
        double sumX1X1 = 0, sumX1X2 = 0, sumX2X2 = 0;
        double sumX1Y = 0, sumX2Y = 0;

        for (int i = 0; i < n; i++) {
            double x1 = data[0][i];
            double x2 = data[1][i];
            double y = data[2][i];

            sumX1 += x1;
            sumX2 += x2;
            sumY += y;

            sumX1X1 += x1 * x1;
            sumX1X2 += x1 * x2;
            sumX2X2 += x2 * x2;
            sumX1Y += x1 * y;
            sumX2Y += x2 * y;
        }

        double[][] matrix = {
                {n, sumX1, sumX2, sumY},
                {sumX1, sumX1X1, sumX1X2, sumX1Y},
                {sumX2, sumX1X2, sumX2X2, sumX2Y}
        };

        GaussMethod(matrix, 3);
        System.out.println("Розвʼязання системи рівнянь:");
        a0 = matrix[0][3];
        a1 = matrix[1][3];
        a2 = matrix[2][3];
        System.out.println("a0 = " + a0);
        System.out.println("a1 = " + a1);
        System.out.println("a2 = " + a2);

        System.out.println("\nЗначення функції у точці (1.5; 3): " + func(1.5, 3));

        double rSquared = R2(data);
        System.out.println("\nКоефіцієнт детермінації R²: " + rSquared);

        plotDataAndEquation(data);

    }

    public static void GaussMethod(double[][] matrix, int size) {
        for (int i = 0; i < size; i++) {
            int maxRow = i;
            for (int k = i + 1; k < size; k++) {
                if (Math.abs(matrix[k][i]) > Math.abs(matrix[maxRow][i])) {
                    maxRow = k;
                }
            }
            double[] temp = matrix[i];
            matrix[i] = matrix[maxRow];
            matrix[maxRow] = temp;

            double diag = matrix[i][i];
            if (diag == 0) {
                System.out.println("Система не має розв’язку або має безліч розв’язків.");
                return;
            }

            for (int k = 0; k < size + 1; k++) {
                matrix[i][k] /= diag;
            }

            for (int j = 0; j < size; j++) {
                if (j == i) continue;
                double factor = matrix[j][i];
                for (int k = 0; k < size + 1; k++) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }
    }

    public static double func (double x1, double x2){
        return (a0 + a1*x1 + a2*x2);
    }

    public static double R2 (double data[][]){
        int n = data[0].length;
        double numerator = 0, denominator = 0;

        double y_mean = Arrays.stream(data[2]).average().orElse(0);

        for (int i = 0; i < n; i++) {
            double x1 = data[0][i];
            double x2 = data[1][i];
            double y = data[2][i];

            numerator += Math.pow((func(x1, x2) - y), 2);
            denominator += Math.pow((y - y_mean), 2);
        }

        return (denominator != 0) ? (1 - numerator / denominator) : Double.NaN;
    }

    public static void plotDataAndEquation(double[][] data) {
        XYSeries series = new XYSeries("Дані");
        for (int i = 0; i < data[0].length; i++) {
            series.add(data[0][i], data[2][i]);
        }

        XYSeries equationSeries = new XYSeries("МНК Рівняння");
        for (double x1 = -1; x1 <= 3; x1 += 0.1) {
            double x2 = 2;
            double y = a0 + a1 * x1 + a2 * x2;
            equationSeries.add(x1, y);
        }

        // Розв'язок рівняння для x1=1.5 та x2=3
        double x1_solution = 1.5;
        double x2_solution = 3;
        double y_solution = a0 + a1 * x1_solution + a2 * x2_solution;

        // Додати точку на графік
        XYSeries solutionPointSeries = new XYSeries("Точка рішення");
        solutionPointSeries.add(x1_solution, y_solution);

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series);
        dataset.addSeries(equationSeries);
        dataset.addSeries(solutionPointSeries);

        JFreeChart chart = ChartFactory.createScatterPlot(
                "Метод найменших квадратів",
                "X1",
                "Y",
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, false);
        renderer.setSeriesShapesVisible(0, true);
        renderer.setSeriesPaint(0, Color.RED);
        renderer.setSeriesLinesVisible(1, true);
        renderer.setSeriesPaint(1, Color.BLUE);
        renderer.setSeriesShapesVisible(2, true); // Для точки рішення
        renderer.setSeriesShape(2, new Ellipse2D.Double(-3, -3, 6, 6)); // Використовуємо коло для точки
        renderer.setSeriesPaint(2, Color.GREEN); // Колір точки рішення

        chart.getXYPlot().setRenderer(renderer);

        JFrame frame = new JFrame("Графік");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

}

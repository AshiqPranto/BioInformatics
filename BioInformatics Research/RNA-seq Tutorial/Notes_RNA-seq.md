### Data quality assessment
Data quality assessment (QA) and exploration are crucial steps in the analysis of gene expression data. These steps are typically performed at the beginning of the analysis process to ensure the reliability and integrity of the data before proceeding to more advanced analyses such as normalization and differential expression testing.

The primary goal of QA is to identify any technical issues or anomalies in the data that could impact the accuracy of downstream analyses. These issues might include sample outliers, batch effects, or other artifacts that can introduce biases and lead to incorrect conclusions.

Here are some common QA steps and techniques, along with an example:

Comparing Replicates: One of the first things to check is the consistency between replicate samples. Replicates are samples that are supposed to be identical, providing a way to assess technical variability. You can visualize the distribution of counts or expression values between replicates to ensure they are similar.
```
boxplot(counts[, c("sample1", "sample2")], main="Replicate Comparison")
```
<B><u>Scatterplot</u></B> of Replicates: Scatterplots are useful for comparing the gene expression levels between replicates. Points should ideally fall along a diagonal line, indicating strong correlation between replicates.
```
plot(counts$sample1, counts$sample2, main="Scatterplot of Replicates")
```
### Skewed distribution

A "skewed distribution" refers to a statistical distribution in which the data points are not evenly distributed around the mean (average) value. Instead, the distribution is stretched or elongated in one direction, often resulting in a tail that extends to one side of the mean more than the other. Skewed distributions can take different shapes:

Positively Skewed (Right Skewed): In a positively skewed distribution, the tail of the distribution extends towards the right side. This means that there are few data points with relatively high values that pull the mean and median in that direction.

Negatively Skewed (Left Skewed): In a negatively skewed distribution, the tail of the distribution extends towards the left side. Here, there are few data points with relatively low values that pull the mean and median in that direction.

In the context of gene expression analysis and RNA-Seq data, a skewed distribution of count values refers to the fact that the majority of genes have low expression levels (low counts), while a relatively small number of genes have very high expression levels (high counts). This is a common characteristic of gene expression data, where most genes are not highly active, but a few may be strongly expressed.

A skewed distribution can pose challenges for statistical analysis and interpretation, as it can affect the assumptions of many statistical tests and models that assume a normal (bell-shaped) distribution. To address this, researchers often apply transformations to the data, such as the log2 transformation, which can help to make the distribution more symmetric and closer to a normal distribution. This transformation can make the data more suitable for various analyses and visualization techniques.
### Why we need a normal distribution?
Having data that approximately follows a normal distribution is often desirable in statistical analysis for several reasons:

1. **Assumption of Many Statistical Tests**: Many commonly used statistical tests and methods, such as t-tests, ANOVA, and linear regression, assume that the data are normally distributed. Violations of this assumption can lead to incorrect results or reduced power of the tests.

2. **Model Assumptions**: In many statistical modeling techniques, the normal distribution is assumed for the residuals (the differences between observed values and predicted values). Deviations from normality can affect the accuracy and reliability of model estimates.

3. **Inference and Confidence Intervals**: Normality is important for making valid statistical inferences, calculating confidence intervals, and conducting hypothesis tests. With a normal distribution, these calculations are well-defined and provide meaningful results.

4. **Statistical Power**: When data follow a normal distribution, statistical tests are more powerful, meaning they are better able to detect true effects when they exist. Non-normality can reduce the power of the tests.

5. **Simplifies Interpretation**: Normal distributions have well-understood and interpretable properties. Deviations from normality can complicate the interpretation of results and make it harder to draw meaningful conclusions.

In the context of gene expression data analysis, applying transformations to the count data, such as the log2 transformation, helps to reduce the skewness and make the data distribution more symmetric. This can make the data more amenable to standard statistical techniques and modeling approaches that assume normality. However, it's important to note that while normality is often preferred, some analyses, like many bioinformatics methods, have been designed to handle count data directly without requiring a strict normal distribution assumption.
### The variance of count data
The variance of count data refers to the spread or dispersion of the counts around their mean value. In other words, it measures how much individual count values deviate from the average count. Variance is a fundamental statistical concept that helps us understand the degree of variability within a set of data points.

Let's consider an example using a hypothetical gene expression dataset, where we have counts of gene expression for three different genes across multiple samples:

| Sample | Gene A | Gene B | Gene C |
|--------|--------|--------|--------|
| Sample 1 | 100 | 200 | 150 |
| Sample 2 | 120 | 180 | 160 |
| Sample 3 | 110 | 210 | 140 |
| Sample 4 | 130 | 190 | 170 |
| Sample 5 | 105 | 220 | 145 |

For each gene, we can calculate the mean count and the variance. Let's calculate the mean and variance of counts for Gene A as an example:

1. Calculate the Mean (Average):
   Mean = (100 + 120 + 110 + 130 + 105) / 5 = 113

2. Calculate the Variance:
   Variance = [(100 - 113)^2 + (120 - 113)^2 + (110 - 113)^2 + (130 - 113)^2 + (105 - 113)^2] / 4
           = [169 + 9 + 9 + 289 + 64] / 4
           = 530 / 4
           = 132.5

The variance of the counts for Gene A is 132.5. This value indicates the average squared difference between each count and the mean count. A higher variance indicates that the counts are more spread out from the mean, while a lower variance indicates that the counts are closer to the mean.

In the context of gene expression analysis, understanding the variance of count data is important because it provides insights into the stability and reliability of the measurements. Genes with high variance may exhibit more dynamic expression patterns, indicating potential biological significance or variability in response to different conditions. On the other hand, low-variance genes may show more consistent expression levels across samples.

Researchers often consider variance when performing quality control, data normalization, and statistical analysis in gene expression studies. Variance stabilization techniques are applied to make data more amenable to statistical tests and to ensure that differences observed are not solely driven by high or low variance in the data.
### The variance of the data increases as the mean increases
The relationship between the variance and the mean of a dataset is a common phenomenon known as the "variance-mean relationship." In many cases, the variance of data tends to increase as the mean of the data increases. This relationship can be observed in various real-world scenarios and has important implications for understanding and analyzing data.

Let's illustrate this relationship with a simple example. Suppose we have two datasets representing the number of daily sales for two different stores over a week:

Store A:
| Day | Sales |
|-----|-------|
| 1 | 100 |
| 2 | 110 |
| 3 | 105 |
| 4 | 120 |
| 5 | 115 |
| 6 | 125 |
| 7 | 130 |

Store B:
| Day | Sales |
|-----|-------|
| 1 | 10 |
| 2 | 15 |
| 3 | 12 |
| 4 | 18 |
| 5 | 14 |
| 6 | 20 |
| 7 | 25 |

Now, let's calculate the mean and variance for each store's sales data:

Store A:
Mean (μA) = (100 + 110 + 105 + 120 + 115 + 125 + 130) / 7 = 114.29
Variance (σ²A) = [(100 - 114.29)² + (110 - 114.29)² + ... + (130 - 114.29)²] / 6 ≈ 133.81

Store B:
Mean (μB) = (10 + 15 + 12 + 18 + 14 + 20 + 25) / 7 = 16.14
Variance (σ²B) = [(10 - 16.14)² + (15 - 16.14)² + ... + (25 - 16.14)²] / 6 ≈ 22.67

In this example, we can observe that the variance of the sales data for Store A (σ²A) is higher than the variance of the sales data for Store B (σ²B), even though the mean sales for Store A (μA) is higher than the mean sales for Store B (μB). This illustrates the variance-mean relationship: as the mean increases, the variance also tends to increase.

The variance-mean relationship is a fundamental concept in statistics and has implications for various data analysis techniques. It's important to be aware of this relationship when interpreting and analyzing datasets, as it can impact the conclusions drawn from statistical tests and modeling.

### BoxPlot
[Video](https://www.youtube.com/watch?v=INSIyaZUXIY)

The boxplot method is a graphical tool used to visualize the distribution of data, in this case, the pseudocounts of gene expression, across different samples. It provides a compact summary of the data's central tendency, spread, and potential outliers.

In a boxplot display:
- A box is drawn that spans from the 25th percentile (lower quartile) to the 75th percentile (upper quartile) of the data distribution. This box represents the interquartile range (IQR), which contains the middle 50% of the data.
- A line (referred to as the "median") is placed within the box at the level that separates the data into two halves: half of the data points are above the median, and half are below.
- Whiskers extend from the edges of the box to the minimum and maximum values within a certain range. The exact length of the whiskers can vary based on statistical considerations.
- Data points that fall beyond the whiskers are typically plotted as individual points, which may indicate potential outliers.

Here's an example of how this would look using pseudocount data for gene expression across different samples:

Suppose we have three samples (S1, S2, S3) and their respective pseudocounts (counts after adding a small constant for stability):
```
Sample   Pseudocounts
S1       10.2
S2       8.5
S3       12.8
```

In a boxplot display of these pseudocounts, we would:
- Draw a box from the 25th percentile (8.5) to the 75th percentile (12.8).
- Place a line (median) at the value 10.2 (the middle value of the sorted data).
- Extend whiskers from the edges of the box to the minimum (8.5) and maximum (12.8) values.
- Plot the data points 8.5 and 12.8 as individual points beyond the whiskers.

The resulting boxplot visually summarizes the distribution of pseudocounts across the three samples, providing insights into the central tendency and variability of the data. It also helps identify potential outliers if any data points fall beyond the whiskers.

The boxplot is a useful tool for quickly assessing the distribution of data and identifying any patterns or discrepancies across different samples.

In this boxplot, the box, median line, and whiskers will be based on the values of the original dataset. However, the outlier value (15.5) will be plotted as an individual point beyond the whiskers. This indicates that the value of 15.5 is an outlier compared to the rest of the data.

Outliers are data points that significantly differ from the majority of the data points. They can sometimes indicate errors, anomalies, or extreme observations in the data. In the context of the boxplot, data points that fall beyond the whiskers are often marked as potential outliers, and they can be further investigated to determine if they are genuine data points or if there are any issues with the measurement or recording process.
### Histogram or density plot
In the context of gene expression analysis, pseudocounts are often used to transform raw count data in order to make it more suitable for various statistical analyses. After applying pseudocounts to the data, it's important to assess the distribution of the transformed values. This can be done using various visualization methods, such as histograms or density plots.

A histogram is a graphical representation of the distribution of data. It divides the range of values into a set of bins and counts how many data points fall into each bin. In the case of pseudocounts, a histogram can provide insights into the frequency of different values after the transformation. Each bin in the histogram represents a range of values, and the height of the bars indicates how many pseudocount values fall within that range.

A density plot is a smoothed version of a histogram. Instead of using discrete bins, a density plot uses a continuous curve to represent the distribution of data. It provides a smoother visualization of the data distribution and can help identify patterns or modes in the distribution.

When assessing pseudocount distributions, histograms or density plots can reveal important characteristics of the transformed data. For example, the presence of multiple modes or peaks in the distribution might indicate underlying subpopulations or patterns in the data. Detecting a secondary mode in the distribution could suggest the presence of a distinct subset of genes with different expression patterns.

In summary, histograms and density plots of pseudocount distributions are valuable tools for exploring and understanding the transformed data, allowing researchers to identify potential features that might not be evident from the raw counts alone.

### PCA
[Tutorial](https://www.youtube.com/watch?v=FgakZw6K1QQ)<br>
Pca creates plots based on correlation among samples<br>
**what is the decision we can take from PCA plot?**
1. Sample Clustering and Separation: PCA can reveal how samples group together based on their overall similarities or differences. Clusters of samples on the PCA plot might indicate underlying patterns, such as different experimental conditions, treatments, or biological groups.
2. Outliers: Outliers or extreme data points can be easily identified in a PCA plot. These points may represent data quality issues, technical artifacts, or true biological variations.
3. Normalization Assessment: PCA can help you assess the effectiveness of normalization methods. If samples from different batches or conditions cluster together on the PCA plot, it might suggest that the normalization was successful in removing batch effects.
**Details explanation of Normalization assesment**
In the context of a PCA plot, normalization assessment involves using the plot to evaluate whether the normalization methods applied to your data were effective in removing batch effects or other sources of technical variability. Here's a broader explanation of how PCA can help with normalization assessment:

1. **Batch Effects:** Batch effects occur when different sets of samples are processed or sequenced at different times or under different conditions. These batch effects can introduce unwanted variability into the data and may obscure true biological differences. In a PCA plot, samples from the same batch or condition should cluster together if batch effects have been successfully removed.

2. **Effectiveness of Normalization:** After applying normalization methods, you can create a PCA plot using the normalized data. If samples from different batches or conditions are no longer clearly separated along the principal components, it suggests that the normalization was effective in mitigating batch effects. Instead, the dominant sources of variation in the data are likely to be biologically meaningful.

3. **Data Reproducibility:** Replicate samples or technical replicates should cluster together closely on the PCA plot, indicating good reproducibility between measurements. If normalization has been successful, replicate samples should be tightly grouped, and their positions on the plot should be relatively consistent.

4. **Interpretation:** If you observe that samples cluster according to their biological conditions rather than by batch or technical factors, it suggests that the normalization has worked well, allowing you to focus on the biological variations of interest.

5. **Next Steps:** If normalization did not effectively remove batch effects, you might need to explore alternative normalization methods or perform additional quality control steps to improve the data quality before downstream analyses.

In summary, PCA can serve as a visual tool to assess the impact of normalization on the distribution of your data and help you determine if normalization has successfully removed unwanted technical variations. Keep in mind that PCA is just one part of the overall quality control process, and additional analyses, statistical tests, and biological validation should be conducted to ensure accurate and reliable results.
### MDS
[Tutorial](https://www.youtube.com/watch?v=GEn-_dAyYME)

MDS creates plots based on distance among samples


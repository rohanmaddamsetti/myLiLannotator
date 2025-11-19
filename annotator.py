import ollama
import numpy as np
import pandas as pd
import re

model_name = "llama3.2:latest"
retry_on = True
small_example = True
n_example = 10

df = pd.read_csv("../data/Maddamsetti2024/gbk-annotation-table.csv", na_filter=False)
ground_truth = pd.read_csv("../data/Maddamsetti2024/computationally-annotated-gbk-annotation-table.csv", na_filter=False)

if small_example:
    df = df.head(n=n_example)
    ground_truth = ground_truth.head(n=n_example)

nrow = df.shape[0]

annotations = np.full(nrow, "NoAnnotation")
categories = ["Agriculture",
"Animals",
"Food",
"Freshwater",
"Human-impacted",
"Humans",
"Livestock",
"Marine",
"Plants",
"Sediment",
"Soil",
"Terrestrial",
"NA"]

categories_str = ', '.join(['"' + category + '"' for category in categories])

for i in range(nrow):
    input_host = df.iloc[i]['host']
    input_isolation_source = df.iloc[i]['isolation_source']

    output = ollama.generate(
        model=model_name,
        prompt='Consider a microbial sample from the host "' + input_host + '" and the isolation source "' + input_isolation_source + '". Please categorize the sample as one of the following categories: ' + categories_str + '. Give a strictly one-word response.'
    )

    if any(re.sub('[."]', '', output['response']) == category for category in categories):
        annotations[i] = re.sub('[."]', '', output['response'])
    elif retry_on:
        output_retry = ollama.chat(
        model=model_name,
        messages=[
            {"role": "system", "content": "You are an annotation tool helping a researcher categorize microbial samples. You will give one-word responses from the following options: " + categories_str + "."},
            {"role": "user", "content": 'Consider a microbial sample from the host "' + input_host + '" and the isolation source "' + input_isolation_source + '". Please categorize the sample as one of the following categories: ' + categories_str + '. Give a strictly one-word response.'},
            {"role": "assistant", "content": output['response']},
            {"role": "user", "content": "That response is not one of the valid options. Please choose the best match from the following options: " + categories_str + '.'}
        ]
        )

        if (any(re.sub('[."]', '', output_retry['message']['content']) == category for category in categories)):
            annotations[i] = re.sub('[."]', '', output_retry['message']['content'])
        else:
            print(output_retry['message']['content'])

    else:
        print(output['response'])



df['Annotation'] = annotations

df.to_csv("../results/latest_result_gbk-annotation-table.csv", index=False)

n_correct = sum(df['Annotation'] == ground_truth['Annotation'])
print("Accuracy: " + str(n_correct/nrow) + " (" + str(n_correct) + " out of " + str(nrow) + ")")

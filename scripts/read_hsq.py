import os
import pandas as pd

def read_hsq(hsq_fn):
    output_records = []
    hsq_file = open(hsq_fn)
    print hsq_fn
    for line in hsq_file:
        line = line.strip().split()
        if line[0] not in ("V(G)/Vp", "V(GxE)/Vp"):
            continue
        else:
            name = 'K_D' if 'GxE' in line[0] else 'K_G'
            output_records.append({"sigma_2": float(line[1]),
                                   "SE": float(line[2]),
                                   "K_type": name})
    return output_records


def read_hsqs(hsq_fns, output_fn):
    output_records = []
    for fn in hsq_fns:
        records = read_hsq(fn)
        for record in records:
            record['phenotype'] = fn.split('___')[-1][:-4]
        output_records += records
    df = pd.DataFrame(output_records)
    df[['SE', 'sigma_2']] *= 100
    df = df.pivot(index='phenotype', columns='K_type')
    # df.columns = df.columns.swaplevel(0,1)
    df.to_csv(output_fn)


if __name__ == "__main__":
    hsq_folder = '../hsq'
    qnormed_hsq_fns = [os.path.join(hsq_folder, x)
                       for x in os.listdir(hsq_folder) if 'qnormed' in x]
    read_hsqs(qnormed_hsq_fns, 'qnormed_variance_components.csv')
    non_qnormed_hsq_fns = [os.path.join(hsq_folder, x)
                           for x in os.listdir(hsq_folder) if 'qnormed' not in x]
    read_hsqs(non_qnormed_hsq_fns, 'non_qnormed_variance_components.csv')


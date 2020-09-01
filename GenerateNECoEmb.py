import argparse
import pandas as pd
import datetime
from gensim.models import Word2Vec
import os


def learn_embeddings(nbhd_file_name, d, w, workers):
    '''
    Learn embeddings by optimizing the Skipgram objective using SGD.
    '''

    walks = read_nbhd(nbhd_file_name)

    model = Word2Vec(walks, size=d, window=w, min_count=0, sg=1, workers=workers,
                     iter=1)
    emb_file_name = os.path.splitext(nbhd_file_name)[0] + ".emb"
    model.wv.save_word2vec_format(emb_file_name)
    print("Finished embedding file: " + emb_file_name)
    return

def read_nbhd(file_name):
    try:
        dfNbhd = pd.read_csv(file_name, sep='|',
                             lineterminator='\n', header=None, dtype=str)
    except os.error as e:
        raise
    dfNbhd = dfNbhd[0].str.split('\t', expand=True)
    nbhd = dfNbhd.astype(str).values.tolist()

    nbhd = [[node for node in walk if node != 'nan' and node != 'None'] for walk in nbhd]
    walks = [list(map(str, walk)) for walk in nbhd]
    return walks

def main(args):
    t0 = datetime.datetime.now().replace(microsecond=0)

    learn_embeddings(args.file, args.dimension, args.window, args.workers)

    print("Time to generate nbhd embeddings: " +
          (str)(datetime.datetime.now().replace(microsecond=0) - t0))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='Neighborhood file name.', default='/users/cdursun/geco/nbhd_ic_Gg_Top250.txt')
    parser.add_argument('-d', '--dimension', help='Latent feature vector size.', type=int, default=512)
    parser.add_argument('-w', '--window', help='Context size for word2vec.', type=int, default=5)
    parser.add_argument('-workers', '--workers', help='Number of processes to be used by word2vec', type=int, default=40)
    args = parser.parse_args()
    main(args)

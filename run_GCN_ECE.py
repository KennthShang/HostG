import  torch
from    torch import nn
from    torch import optim
from    torch.nn import functional as F

import  numpy as np
from    data import load_data, preprocess_features, preprocess_adj, sample_mask
import  model
from    config import  args
from    utils import masked_loss, masked_acc, masked_ECE
import  pickle as pkl
import  scipy.sparse as sp
import argparse
from scipy.special import softmax
from sklearn.metrics import classification_report
from collections import Counter

import random


################################################################################
###############################  Input Params  #################################
################################################################################

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--gpus', type=int, default = 0)
parser.add_argument('--t', type=float, default=0.0)
inputs = parser.parse_args()


seed = 123
np.random.seed(seed)
torch.random.manual_seed(seed)


if torch.cuda.is_available():
    torch.cuda.set_device(inputs.gpus)
else:
    print("Running with cpu")

################################################################################
############################  Loading dataset  #################################
################################################################################

def GCN(taxa):
    adj        = pkl.load(open("GCN_data/"+taxa+"_contig.graph",'rb'))
    labels     = pkl.load(open("GCN_data/"+taxa+"_contig.label",'rb'))
    features   = pkl.load(open("GCN_data/"+taxa+"_contig.feature",'rb'))
    idx_test   = pkl.load(open("GCN_data/"+taxa+"_contig.test_id",'rb'))


    idx_test = np.array(idx_test)
    labels = np.array(labels)
    y_train = np.zeros(labels.shape)
    y_test = np.zeros(labels.shape)


    idx_train = np.array([i for i in range(len(labels))])
    idx_train = np.array([i for i in range(len(labels)) if i not in idx_test])


    train_mask = sample_mask(idx_train, labels.shape[0])
    test_mask = sample_mask(idx_test, labels.shape[0])


    y_train[train_mask] = labels[train_mask]
    y_test[test_mask] = labels[test_mask]


    features = sp.csc_matrix(features)

    print('adj:', adj.shape)
    print('features:', features.shape)
    print('y:', y_train.shape, y_test.shape) # y_val.shape, 
    print('mask:', train_mask.shape, test_mask.shape) # val_mask.shape


    ################################################################################
    ###########################  Data preprocessing  ###############################
    ################################################################################

    features = preprocess_features(features) # [49216, 2], [49216], [2708, 1433]
    supports = preprocess_adj(adj)

    if torch.cuda.is_available():
        device = torch.device('cuda')
        train_label = torch.from_numpy(y_train).long().to(device)
        num_classes = max(labels)+1
        train_mask = torch.from_numpy(train_mask.astype(bool)).to(device)
        test_label = torch.from_numpy(y_test).long().to(device)
        test_mask = torch.from_numpy(test_mask.astype(bool)).to(device)
        # graph
        i = torch.from_numpy(features[0]).long().to(device)
        v = torch.from_numpy(features[1]).to(device)
        feature = torch.sparse.FloatTensor(i.t(), v, features[2]).float().to(device)
        feature = feature.to_dense()
        i = torch.from_numpy(supports[0]).long().to(device)
        v = torch.from_numpy(supports[1]).to(device)
        support = torch.sparse.FloatTensor(i.t(), v, supports[2]).float().to(device)
        support = support.to_dense()
    else:
        train_label = torch.from_numpy(y_train).long()
        num_classes = max(labels)+1
        train_mask = torch.from_numpy(train_mask.astype(bool))
        test_label = torch.from_numpy(y_test).long()
        test_mask = torch.from_numpy(test_mask.astype(bool))
        # graph
        i = torch.from_numpy(features[0]).long()
        v = torch.from_numpy(features[1])
        feature = torch.sparse.FloatTensor(i.t(), v, features[2]).float()
        feature = feature.to_dense()
        i = torch.from_numpy(supports[0]).long()
        v = torch.from_numpy(supports[1])
        support = torch.sparse.FloatTensor(i.t(), v, supports[2]).float()
        support = support.to_dense()


    print('x :', feature)
    print('sp:', support)
    #num_features_nonzero = feature._nnz()
    feat_dim = feature.shape[1]


    ################################################################################
    ##############################  Training model  ################################
    ################################################################################

    def accuracy(out, mask):
        pred = np.argmax(out, axis = 1)
        mask_pred = np.array([pred[i] for i in range(len(labels)) if mask[i] == True])
        mask_label = np.array([labels[i] for i in range(len(labels)) if mask[i] == True])
        return np.sum(mask_label == mask_pred)/len(mask_pred)


    node_dim = adj.shape[0]
    net = model.GCN(node_dim, feat_dim, num_classes, 0)
    if torch.cuda.is_available():
        net.to(device)

    optimizer = optim.Adam(net.parameters(), lr=0.01)#args.learning_rate



    _ = net.train()
    for epoch in range(args.epochs*2):
        # forward pass
        out = net((feature, support))
        loss = masked_loss(out, train_label, train_mask)
        #loss += masked_ECE(out, train_label, train_mask)
        loss += args.weight_decay * net.l2_loss()
        # backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        # output
        if epoch % 10 == 0:
            # calculating the acc
            _ = net.eval()
            out = net((feature, support))
            if torch.cuda.is_available():
                acc_train = accuracy(out.detach().cpu().numpy(), train_mask.detach().cpu().numpy())
            else:
                acc_train = accuracy(out.detach().numpy(), train_mask.detach().numpy())
            print("Traing epoch: {0} || training loss: {1} || training Acc {2}".format(epoch, round(loss.item(), 2), round(acc_train, 2)))
        _ = net.train()





    ################################################################################
    #################################  Prediction   ################################
    ################################################################################

    net.eval()
    out = net((feature, support))
    if torch.cuda.is_available():
        out = out.cpu().detach().numpy()
    else:
        out = out.detach().numpy()

    pred = np.argmax(out, axis = 1)


    label2int   = pkl.load(open("GCN_data/"+taxa+"_label2int.dict",'rb'))
    id2node     = pkl.load(open("GCN_data/"+taxa+"_id2node.dict",'rb'))

    int2label ={idx: label for label, idx in label2int.items()}

    with open("tmp_pred/prediction_"+taxa+".csv", 'w') as f_out:
        _ = f_out.write("contig_names,"+ taxa +"\n")
        for idx, node in id2node.items():
            if "cherry" in node:
                if labels[idx] != -1:
                    _ = f_out.write(node + "," + int2label[labels[idx]] + "\n")
                    #print(node + "," + int2label[labels[idx]])
                else:
                    if max(softmax(out[idx])) > inputs.t:
                        _ = f_out.write(node + "," + int2label[pred[idx]] + "\n")
                        #print(node + "," + int2label[pred[idx]])


taxa_list = ["phylum", 'class', 'order', 'family', 'genus']
for taxa in taxa_list:
    GCN(taxa)

import  torch
from    torch import nn
from    torch.nn import functional as F


def masked_ECE(out, train_label, train_mask):
    device = torch.device('cuda')
    m = nn.Softmax(dim=1)
    prob_box = dict.fromkeys([0, 1, 2, 3, 4])
    pred_box = dict.fromkeys([0, 1, 2, 3, 4])
    label_box = dict.fromkeys([0, 1, 2, 3, 4])
    # calculate pred
    prob = torch.max(m(out), 1)[0]
    pred = torch.argmax(out, 1)
    # calculate box
    for i in range(len(prob)):
        if train_mask[i] == True:
            if prob[i] > 0.8:
                try:
                    prob_box[4].append(prob[i])
                    pred_box[4].append(pred[i])
                    label_box[4].append(train_label[i])
                except:
                    prob_box[4] = [prob[i]]
                    pred_box[4] = [pred[i]]
                    label_box[4]= [train_label[i]]
            elif prob[i] > 0.6:
                try:
                    prob_box[3].append(prob[i])
                    pred_box[3].append(pred[i])
                    label_box[3].append(train_label[i])
                except:
                    prob_box[3] = [prob[i]]
                    pred_box[3] = [pred[i]]
                    label_box[3]= [train_label[i]]
            elif prob[i] > 0.4:
                try:
                    prob_box[2].append(prob[i])
                    pred_box[2].append(pred[i])
                    label_box[2].append(train_label[i])
                except:
                    prob_box[2] = [prob[i]]
                    pred_box[2] = [pred[i]]
                    label_box[2]= [train_label[i]]
            elif prob[i] > 0.2:
                try:
                    prob_box[1].append(prob[i])
                    pred_box[1].append(pred[i])
                    label_box[1].append(train_label[i])
                except:
                    prob_box[1] = [prob[i]]
                    pred_box[1] = [pred[i]]
                    label_box[1]= [train_label[i]]
            else:
                try:
                    prob_box[0].append(prob[i])
                    pred_box[0].append(pred[i])
                    label_box[0].append(train_label[i])
                except:
                    prob_box[0] = [prob[i]]
                    pred_box[0] = [pred[i]]
                    label_box[0]= [train_label[i]]
    # calculate ECE
    ECE = 0
    for key in label_box.keys():
        if label_box[key] == None:
            continue
        accuarcy = torch.sum(torch.FloatTensor(label_box[key]) == torch.FloatTensor(pred_box[key]))/len(pred_box[key])
        confidence = torch.mean(torch.FloatTensor(prob_box[key]))
        gap = len(pred_box[key])*torch.abs(accuarcy-confidence)
        ECE+=gap
    return ECE.to(device)/len(out)


def masked_loss(out, label, mask):
    #if torch.cuda.is_available():
    #    w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0]).cuda()
    #else:
    #    w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0])
    #loss = F.cross_entropy(out, label, w, reduction='none')
    loss = F.cross_entropy(out, label, reduction='none')
    #all phage
    #w = torch.Tensor([3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 3.0, 2.0, 3.0]).cuda()
    #loss = F.cross_entropy(out, label, w, reduction='none')
    mask = mask.float()
    mask = mask / mask.mean()
    loss *= mask
    loss = loss.mean()
    return loss


def masked_acc(out, label, mask):
    # [node, f]
    pred = out.argmax(dim=1)
    correct = torch.eq(pred, label).float()
    mask = mask.float()
    mask = mask / mask.mean()
    correct *= mask
    acc = correct.mean()
    return acc



def sparse_dropout(x, rate, noise_shape):
    """

    :param x:
    :param rate:
    :param noise_shape: int scalar
    :return:
    """
    random_tensor = 1 - rate
    random_tensor += torch.rand(noise_shape).to(x.device)
    dropout_mask = torch.floor(random_tensor).byte()
    i = x._indices() # [2, 49216]
    v = x._values() # [49216]

    # [2, 4926] => [49216, 2] => [remained node, 2] => [2, remained node]
    i = i[:, dropout_mask]
    v = v[dropout_mask]

    out = torch.sparse.FloatTensor(i, v, x.shape).to(x.device)

    out = out * (1./ (1-rate))

    return out


def dot(x, y, sparse=False):
    if sparse:
        res = torch.sparse.mm(x, y)
    else:
        res = torch.mm(x, y)

    return res


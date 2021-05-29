m = nn.Softmax(dim=1)
# calculate pred
prob = torch.max(m(out), 1)[0]
pred = torch.argmax(out, 1)
test_label

test_prob = []
test_acc = []
prob2acc = {}


for i in range(len(prob)):
    if test_mask[i] == True:
        tmp_prob = prob[i].cpu().detach().numpy().tolist()
        tmp_pred = pred[i].cpu().detach().numpy().tolist()
        tmp_label= test_label[i].cpu().detach().numpy().tolist()
        try:
            prob2acc[tmp_prob].append(tmp_label == tmp_pred)
        except:
            prob2acc[tmp_prob] = [tmp_label == tmp_pred]


keys = prob2acc.keys()
sorted_prob2acc = {sorted_key: prob2acc[sorted_key] for sorted_key in sorted(keys, reverse=True)}

cnt = 0
correct = 0
for key in sorted_prob2acc:
    tmp_pred_list = sorted_prob2acc[key]
    cnt += len(tmp_pred_list)
    correct += np.sum(tmp_pred_list)
    print("recall rate: {0} || accuarcy: {1}".format(cnt/70, correct/cnt))
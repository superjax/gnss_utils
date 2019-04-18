import numpy as np
import matplotlib.pyplot as plt

DataType = np.dtype([
	("t", np.float64),
	("lla", (np.float64, 3)),
	("ref_lla", (np.float64, 3)),
	("vel", (np.float64, 3)),
	("ref_vel", (np.float64, 3))
	])

x = np.fromfile("/tmp/positioning_demo.bin", dtype=DataType)

plt.figure()
plt.suptitle("lla")
for i in range(3):
	plt.subplot(3,1,i+1)
	plt.plot(x["t"], x["lla"][:,i], label="lla")
	plt.plot(x["t"], x["ref_lla"][:,i], label="ref")
	if i == 0:
		plt.legend()

plt.figure()
plt.suptitle("vel")
for i in range(3):
	plt.subplot(3,1,i+1)
	plt.plot(x["t"], x["vel"][:,0], label="vel")
	plt.plot(x["t"], x["ref_vel"][:,0], label="ref")
	if i == 0:
		plt.legend()


plt.show()

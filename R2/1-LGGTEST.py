# -*- coding: utf-8 -*-
import os
import sys
import random
import math
import numpy as np
import skimage.io
import matplotlib
import matplotlib.pyplot as plt
import cv2
import time
from mrcnn.config import Config
from datetime import datetime
import radiomics
import radiomics.featureextractor as FEE
import json
# Root directory of the project
ROOT_DIR = os.path.abspath("G:\\LGG/")

# Import Mask RCNN
sys.path.append(ROOT_DIR)  # To find local version of the library
from mrcnn import utils
import mrcnn.model as modellib
from mrcnn import visualize

# Import COCO config
# sys.path.append(os.path.join(ROOT_DIR, "samples/coco/"))  # To find local version
# from samples.coco import coco


# Directory to save logs and trained model
MODEL_DIR = os.path.join(ROOT_DIR, "logs")

# Local path to trained weights file
COCO_MODEL_PATH = os.path.join(MODEL_DIR, "mask_rcnn_shapes_0010.h5")
# Download COCO trained weights from Releases if needed
if not os.path.exists(COCO_MODEL_PATH):
    utils.download_trained_weights(COCO_MODEL_PATH)
    print("cuiwei***********************")

# Directory of images to run detection on
IMAGE_DIR = os.path.join(ROOT_DIR, "test")


class ShapesConfig(Config):
    """Configuration for training on the toy shapes dataset.
    Derives from the base Config class and overrides values specific
    to the toy shapes dataset.
    """
    # Give the configuration a recognizable name
    NAME = "shapes"

    # Train on 1 GPU and 8 images per GPU. We can put multiple images on each
    # GPU because the images are small. Batch size is 8 (GPUs * images/GPU).
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1

    # Number of classes (including background)
    NUM_CLASSES = 1 + 1  # background + 1 shapes

    # Use small images for faster training. Set the limits of the small side
    # the large side, and that determines the image shape.
    IMAGE_MIN_DIM = 320
    IMAGE_MAX_DIM = 384

    # Use smaller anchors because our image and objects are small
    RPN_ANCHOR_SCALES = (8 * 6, 16 * 6, 32 * 6, 64 * 6, 128 * 6)  # anchor side in pixels

    # Reduce training ROIs per image because the images are small and have
    # few objects. Aim to allow ROI sampling to pick 33% positive ROIs.
    TRAIN_ROIS_PER_IMAGE = 100

    # Use a small epoch since the data is simple
    STEPS_PER_EPOCH = 100

    # use small validation steps since the epoch is small
    VALIDATION_STEPS = 50


# import train_tongue
# class InferenceConfig(coco.CocoConfig):
class InferenceConfig(ShapesConfig):
    # Set batch size to 1 since we'll be running inference on
    # one image at a time. Batch size = GPU_COUNT * IMAGES_PER_GPU
    GPU_COUNT = 1
    IMAGES_PER_GPU = 1


config = InferenceConfig()

#model = modellib.MaskRCNN(mode="inference", model_dir=MODEL_DIR, config=config)

# Create model object in inference mode.
model = modellib.MaskRCNN(mode="inference", model_dir=MODEL_DIR, config=config)

# Load weights trained on MS-COCO
model.load_weights(COCO_MODEL_PATH, by_name=True)

# COCO Class names
# Index of the class in the list is its ID. For example, to get ID of
# the teddy bear class, use: class_names.index('teddy bear')
class_names = ['backgroud', 'cancer']
# Load a random image from the images folder
samples_names = next(os.walk(IMAGE_DIR))[1]
for sample_name in samples_names:
    sample_name_path = os.path.join(IMAGE_DIR, sample_name)
    file_names = next(os.walk(sample_name_path))[2]
    for file_name in file_names:
        file_names_path = os.path.join(IMAGE_DIR, sample_name, file_name)
        file_names_outpath_1 = os.path.join(ROOT_DIR, "mask", sample_name)
        if not os.path.exists(file_names_outpath_1):
            os.mkdir(file_names_outpath_1)
        file_names_outpath = os.path.join(file_names_outpath_1, file_name)
        file_names_outpath_masked_1 = os.path.join(ROOT_DIR, "masked", sample_name)
        if not os.path.exists(file_names_outpath_masked_1):
            os.mkdir(file_names_outpath_masked_1)
        file_names_outpath_masked = os.path.join(file_names_outpath_masked_1, file_name)
        file_names_outpath_feature_1 = os.path.join(ROOT_DIR, "features", sample_name)
        if not os.path.exists(file_names_outpath_feature_1):
            os.mkdir(file_names_outpath_feature_1)

        image = cv2.imread(file_names_path)
        a = datetime.now()
    # Run detection
        results = model.detect([image], verbose=1)
        b = datetime.now()
    # Visualize results
        print("shijian", (b - a).seconds)
        r = results[0]
    visualize.display_instances(image, r['rois'], r['masks'], r['class_ids'],
                                    class_names, r['scores'])


        #if len(r['class_ids'])==1 and r['scores'].item(0)>=0.97 and ((r['rois'].item(2) - r['rois'].item(0))*(r['rois'].item(3) - r['rois'].item(1))) < 3000:
        if len(r['class_ids']) == 1:
            #cv2.imshow('img',image)
            n3 = r['rois'].item(3)
            n2 = r['rois'].item(2)
            n1 = r['rois'].item(1)
            n0 = r['rois'].item(0)
            #if r['rois'].item(0) < r['rois'].item(1) and r['rois'].item(2) < r['rois'].item(3):
            #    d = np.array([[[n0 + 10, n3 - 10], [n1 + 30, n3 - 10], [n1 + 30, n2 - 30], [n0 + 10, n2 - 30]]],
            #                 dtype=np.int32)
            #elif r['rois'].item(0) < r['rois'].item(1) and r['rois'].item(2) > r['rois'].item(3):
            #   d = np.array([[[n1, n0], [n3, n0], [n3, n2], [n1, n2]]], dtype=np.int32)
            d = np.array([[[n1, n0], [n3, n0], [n3, n2], [n1, n2]]], dtype=np.int32)
            im = np.zeros(image.shape[:2], dtype="uint8")
            cv2.polylines(im, d, 1, 255)
            cv2.fillPoly(im, d, 255)
            # cv2.imshow("image", image)
            mask = im
            # cv2.imshow("Mask", mask)
            masked = cv2.bitwise_and(image, image, mask=mask)
            # cv2.imshow("Mask to Image", masked)
            # cv2.waitKey(0)
            cv2.imwrite(file_names_outpath, mask, [int(cv2.IMWRITE_JPEG_QUALITY), 100])
            cv2.imwrite(file_names_outpath_masked, masked, [int(cv2.IMWRITE_JPEG_QUALITY), 100])

            """
            Pyradiomics 
            """
                # 文件全部路径
            ori_path = file_names_path
            lab_path = file_names_outpath
            para_path = r'E:\Mask_RCNN_master\pyradiomics\Params.yaml'
            print("originl path: " + ori_path)
            print("label path: " + lab_path)
            print("parameter path: " + para_path)

                # 使用配置文件初始化特征抽取器
            extractor = FEE.RadiomicsFeatureExtractor(para_path)
            print("Extraction parameters:\n\t", extractor.settings)
            print("Enabled filters:\n\t", extractor.enabledImagetypes)
            print("Enabled features:\n\t", extractor.enabledFeatures)

                # 运行
            result = extractor.execute(ori_path, lab_path)  # 抽取特征
            print("Result type:", type(result))  # result is returned in a Python ordered dictionary
            print("")
            print("Calculated features")
            for key, value in result.items():  # 输出特征
                feature_file_name = file_name.split(".")[0] + ".txt"
                f = os.path.join(file_names_outpath_feature_1, feature_file_name)
                teb = len(result)
                with open(f, "a") as file:  # 只需要将之前的”w"改为“a"即可，代表追加内容
                    file.write(str(key) + "\t" + str(value) + "\n")

                # print("\t", key, ":", value)